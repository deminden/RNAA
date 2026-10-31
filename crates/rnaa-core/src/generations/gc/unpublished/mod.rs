//! Reversible cleanup of tracked, terminal unpublished attempts; never result adoption.
mod catalog;
mod files;
mod protection;
#[cfg(test)]
mod tests;
mod tree;
use super::*;
use catalog::{Descriptor, GenerationRow, Receipt};
use protection::Protection;

#[derive(Serialize)]
pub struct UnpublishedCleanupEntry {
    pub generation_id: String,
    pub path: PathBuf,
    pub state: String,
    pub action: String,
    pub reasons: Vec<String>,
    pub logical_bytes: Option<u64>,
}
#[derive(Serialize)]
pub struct UnpublishedCleanupReport {
    pub schema_version: u32,
    pub applied: bool,
    pub busy: bool,
    pub entries: Vec<UnpublishedCleanupEntry>,
    pub errors: Vec<String>,
}
impl UnpublishedCleanupReport {
    fn new(apply: bool) -> Self {
        Self {
            schema_version: 1,
            applied: apply,
            busy: false,
            entries: Vec::new(),
            errors: Vec::new(),
        }
    }
}
impl Database {
    pub fn unpublished_generation_status(&self) -> Result<serde_json::Value> {
        let _pin = super::super::pins::read_lock(self.path())?;
        let rows=Receipt::all(self)?.into_values().map(|row|serde_json::json!({"generation_id":row.descriptor.generation.id,"state":row.state,"operation_id":row.descriptor.operation_id,"trashed_at_ms":row.descriptor.trashed_at_ms,"purge_after_ms":row.descriptor.purge_after_ms})).collect::<Vec<_>>();
        Ok(serde_json::json!({"schema_version":1,"trash":rows}))
    }
    pub fn unpublished_generation_gc(
        &self,
        options: GenerationGcOptions,
    ) -> Result<UnpublishedCleanupReport> {
        let mut report = UnpublishedCleanupReport::new(options.apply);
        let Some(_lock) = super::super::pins::maintenance_lock(self)? else {
            report.busy = true;
            return Ok(report);
        };
        let mut db = self.clone();
        db.generation_maintenance = true;
        let before = Receipt::all(&db)?;
        // Only pre-existing trash may be purged by this invocation. Recovery of
        // an older intent never enrolls freshly discovered paths into that intent.
        if options.apply {
            db.recover_unpublished(&mut report.errors)?;
        }
        let protection = match Protection::capture(&db) {
            Ok(value) => Some(value),
            Err(error) => {
                report.errors.push(format!("protection catalog: {error:#}"));
                None
            }
        };
        let receipts = Receipt::all(&db)?;
        let generations = GenerationRow::unpublished(&db)?;
        let mut plans = Vec::new();
        let mut purges = Vec::new();
        let now = now_ms()?;
        let deadline = now
            .checked_add(i64::try_from(
                options
                    .trash_days
                    .checked_mul(86_400_000)
                    .context("unpublished trash retention overflow")?,
            )?)
            .context("unpublished trash deadline overflow")?;
        for row in generations.values().filter(|row| row.state != "COMMITTED") {
            let mut entry = UnpublishedCleanupEntry {
                generation_id: row.id.clone(),
                path: row.payload.clone(),
                state: row.state.clone(),
                action: String::new(),
                reasons: Vec::new(),
                logical_bytes: None,
            };
            let result = (|| -> Result<()> {
                if let Some(receipt) = receipts.get(&row.id) {
                    entry.state = receipt.state.clone();
                    entry.logical_bytes =
                        Some(receipt.descriptor.tree.iter().map(|item| item.bytes).sum());
                    entry.action = match receipt.state.as_str() {
                        "PURGED" => "already_purged",
                        "TRASHED" => {
                            entry.reasons = if let Some(protection) = &protection {
                                protection.reasons(&db, &row.id)?
                            } else {
                                vec!["unverified_protection_catalog".into()]
                            };
                            if options.purge_expired
                                && before.contains_key(&row.id)
                                && receipt.descriptor.purge_after_ms <= now
                                && entry.reasons.is_empty()
                            {
                                tree::verify(
                                    &tree::container(&db, &row.id, true, false)?,
                                    &receipt.descriptor.tree,
                                    false,
                                )?;
                                purges.push(row.id.clone());
                                "purge"
                            } else {
                                "held_trash_window_or_reference"
                            }
                        }
                        _ => "resume_pending",
                    }
                    .into();
                    return Ok(());
                }
                row.verify_unpublished(&db)?;
                entry.reasons = if let Some(protection) = &protection {
                    protection.reasons(&db, &row.id)?
                } else {
                    vec!["unverified_protection_catalog".into()]
                };
                if !entry.reasons.is_empty() {
                    entry.action = "held_referenced".into();
                    return Ok(());
                }
                let source = tree::container(&db, &row.id, false, false)?;
                if tree::container(&db, &row.id, true, false)?.exists() {
                    bail!("unpublished trash destination already exists");
                }
                let observed = tree::observe(&source, false)?;
                entry.logical_bytes = Some(observed.iter().map(|item| item.bytes).sum());
                let mut catalog = row.clone();
                catalog.state = "ABANDONED".into();
                plans.push(Descriptor {
                    schema_version: 1,
                    operation_id: uuid::Uuid::new_v4().to_string(),
                    generation: catalog,
                    original_state: row.state.clone(),
                    tree: observed,
                    trashed_at_ms: now,
                    purge_after_ms: deadline,
                });
                entry.action = "trash".into();
                Ok(())
            })();
            if let Err(error) = result {
                entry.action = "held_invalid".into();
                entry.reasons.push(format!("{error:#}"));
                report.errors.push(format!("{}: {error:#}", row.id));
            }
            report.entries.push(entry);
        }
        db.report_unindexed(&mut report)?;
        if options.apply {
            db.mark_unpublished(&plans)?;
            for id in purges {
                receipts[&id].transition(&db, "PURGING")?;
            }
            db.recover_unpublished(&mut report.errors)?;
            let after = Receipt::all(&db)?;
            for entry in &mut report.entries {
                if let Some(row) = after.get(&entry.generation_id) {
                    entry.state = row.state.clone();
                }
            }
        }
        Ok(report)
    }
    pub fn restore_unpublished_generation(&self, id: &str) -> Result<()> {
        let _lock = super::super::pins::maintenance_lock(self)?
            .context("project is busy; retry unpublished restore")?;
        let mut db = self.clone();
        db.generation_maintenance = true;
        let pending_restore = Receipt::all(&db)?
            .get(id)
            .is_some_and(|row| row.state == "RESTORING");
        let mut errors = Vec::new();
        db.recover_unpublished(&mut errors)?;
        let rows = Receipt::all(&db)?;
        if pending_restore && !rows.contains_key(id) {
            return Ok(());
        }
        let row = rows
            .get(id)
            .context("no unpublished trash receipt for this generation")?;
        if row.state != "TRASHED" {
            bail!("unpublished generation is not restorable: {}", row.state);
        }
        tree::verify(
            &tree::container(&db, id, true, false)?,
            &row.descriptor.tree,
            false,
        )?;
        row.transition(&db, "RESTORING")?;
        db.recover_unpublished(&mut errors)?;
        if Receipt::all(&db)?.contains_key(id) {
            bail!("unpublished restore is pending: {}", errors.join("; "));
        }
        Ok(())
    }
    fn report_unindexed(&self, report: &mut UnpublishedCleanupReport) -> Result<()> {
        let root = tree::checked_directory(self, &[".rnaa_generations"], false)?;
        if !root.exists() {
            return Ok(());
        }
        let indexed = self.read(|conn| {
            let mut query =
                conn.prepare("SELECT generation_id FROM artifact_generations WHERE project_id=?1")?;
            Ok(query
                .query_map([self.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<BTreeSet<_>, _>>()?)
        })?;
        for entry in fs::read_dir(root)? {
            let entry = entry?;
            let id = entry.file_name().to_string_lossy().into_owned();
            if !indexed.contains(&id) {
                report.entries.push(UnpublishedCleanupEntry {
                    generation_id: id,
                    path: entry.path(),
                    state: "UNINDEXED".into(),
                    action: "held_unindexed".into(),
                    reasons: vec!["no_published_attempt_ownership_witness".into()],
                    logical_bytes: None,
                });
            }
        }
        report
            .entries
            .sort_by(|a, b| a.generation_id.cmp(&b.generation_id));
        Ok(())
    }
}
