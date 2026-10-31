use super::*;
use std::fs::File;
impl Database {
    pub(super) fn mark_unpublished(&self, plans: &[Descriptor]) -> Result<()> {
        for plan in plans {
            let source = tree::container(self, &plan.generation.id, false, false)?;
            if tree::observe(&source, true)? != plan.tree {
                bail!("unpublished files changed before trash marking");
            }
            let path = catalog::witness_path(self, &plan.operation_id, true)?;
            crate::util::publish_file(&path, |temporary| {
                fs::write(temporary, serde_json::to_vec(plan)?)?;
                Ok(())
            })?;
        }
        self.atomic(|db| {
            let rows=GenerationRow::unpublished(db)?;
            for plan in plans {
                let mut expected=plan.generation.clone();expected.state=plan.original_state.clone();
                if rows.get(&expected.id)!=Some(&expected) {bail!("unpublished catalog changed before trash marking");}
                expected.verify_unpublished(db)?;
                let json=serde_json::to_string(plan)?;
                db.write(|conn| {
                    conn.execute("UPDATE artifact_generations SET state='ABANDONED' WHERE generation_id=?1",[&expected.id])?;
                    conn.execute("INSERT INTO unpublished_generation_trash(generation_id,project_id,state,descriptor_json,descriptor_sha256,trashed_at_ms,purge_after_ms) VALUES (?1,?2,'TRASHING',?3,?4,?5,?6)",params![expected.id,db.project_id(),json,digest(&json),plan.trashed_at_ms,plan.purge_after_ms])?;
                    Ok(())
                })?;
                db.append_event("cleanup",None,"unpublished generation trash marked",serde_json::json!({"generation_id":expected.id,"attempt_id":expected.attempt,"operation_id":plan.operation_id,"purge_after_ms":plan.purge_after_ms}))?;
            }
            Ok(())
        })
    }
    pub(super) fn recover_unpublished(&self, errors: &mut Vec<String>) -> Result<()> {
        let rows = Receipt::all(self)?;
        // Restore first: moving a diagnostic container back may resolve aliases.
        // Then build one protected catalog view for the remaining moves/purges,
        // rather than rescanning the complete generation graph per candidate.
        for row in rows.values().filter(|row| row.state == "RESTORING") {
            if let Err(error) = self.recover_unpublished_row(row, None) {
                errors.push(format!(
                    "{} pending restore: {error:#}",
                    row.descriptor.generation.id
                ));
            }
        }
        let pending = rows
            .values()
            .filter(|row| matches!(row.state.as_str(), "TRASHING" | "PURGING"))
            .collect::<Vec<_>>();
        if pending.is_empty() {
            return Ok(());
        }
        let protection = Protection::capture(self);
        for row in pending {
            let result = match &protection {
                Ok(protection) => self.recover_unpublished_row(row, Some(protection)),
                Err(error) => Err(anyhow::anyhow!("protection catalog: {error:#}")),
            };
            if let Err(error) = result {
                errors.push(format!(
                    "{} pending {}: {error:#}",
                    row.descriptor.generation.id, row.state
                ));
            }
        }
        Ok(())
    }
    fn recover_unpublished_row(
        &self,
        row: &Receipt,
        protection: Option<&Protection>,
    ) -> Result<()> {
        let id = &row.descriptor.generation.id;
        let restore = row.state == "RESTORING";
        if !restore {
            let reasons = protection
                .context("unpublished recovery requires a protection snapshot")?
                .reasons(self, id)?;
            if !reasons.is_empty() {
                bail!(
                    "unpublished data acquired references: {}",
                    reasons.join(", ")
                );
            }
        }
        if row.state == "PURGING" {
            let root = tree::container(self, id, true, false)?;
            if root.exists() {
                tree::verify(&root, &row.descriptor.tree, true)?;
                fs::remove_dir_all(&root)?;
                File::open(root.parent().unwrap())?.sync_all()?;
            }
            row.transition(self, "PURGED")?;
        } else if matches!(row.state.as_str(), "TRASHING" | "RESTORING") {
            let source = tree::container(self, id, restore, false)?;
            let target = tree::container(self, id, !restore, true)?;
            match (source.exists(), target.exists()) {
                (true, false) => {
                    tree::verify(&source, &row.descriptor.tree, false)?;
                    tree::sync_move(&source, &target)?;
                }
                (false, true) => tree::verify(&target, &row.descriptor.tree, false)?,
                _ => bail!("ambiguous unpublished cleanup/restore paths"),
            }
            if restore {
                self.atomic(|db| {
                    db.write(|conn| {conn.execute("DELETE FROM unpublished_generation_trash WHERE generation_id=?1 AND state='RESTORING'",[id])?;Ok(())})?;
                    db.append_event("cleanup",None,"unpublished diagnostic data restored without publishing a result",serde_json::json!({"generation_id":id}))
                })?;
            } else {
                row.transition(self, "TRASHED")?;
            }
        } else {
            bail!("unsupported unpublished cleanup state");
        }
        Ok(())
    }
}
