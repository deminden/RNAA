use super::*;
use std::fs::File;

impl Database {
    pub(super) fn mark_file_trash(&self, values: &[Descriptor]) -> Result<()> {
        for value in values {
            self.verify_trash_file(&self.original_file_path(value, false)?, value)?;
            if self.copy_staging_path(value, false)?.exists() {
                bail!("cleanup copy staging identity already exists");
            }
            let path = self.file_receipt_path(&value.trash_id, true)?;
            if path.exists() {
                bail!("file cleanup receipt identity already exists");
            }
            write_json_pretty(&path, value)?;
        }
        self.atomic(|db| {
            for value in values {
                let current=db.list_artifacts_for_run(Some(&value.run))?.into_iter().find(|record|record.path==value.record.path).context("cleanup registration disappeared")?;
                if serde_json::to_value(current)?!=serde_json::to_value(&value.record)? || db.generation_head(OutputGroup::Quant,&value.run)?.is_none_or(|head|head.generation_id!=value.consumer) {
                    bail!("cleanup count or input registration changed before publication");
                }
                if let Some(binding)=&value.binding && db.generation_input_binding(OutputGroup::Quant,&value.run)?.is_none_or(|current|current.binding_id!=*binding) {bail!("cleanup input binding changed");}
                let json=serde_json::to_string(value)?;
                db.write(|conn| {
                    conn.execute("INSERT INTO file_trash(trash_id,project_id,run_accession,canonical_source,operation_id,state,descriptor_json,descriptor_sha256,created_at_ms,purge_after_ms) VALUES (?1,?2,?3,?4,?5,'MOVING',?6,?7,?8,?9)",
                        params![value.trash_id,db.project_id(),value.run,value.canonical_source,value.operation_id,json,digest(json.as_bytes()),value.created_at_ms,value.purge_after_ms])?;
                    Ok(())
                })?;
                db.delete_artifact(&value.record.path)?;
                db.append_event("cleanup",Some(&value.run),"input file trash marked",serde_json::json!({"trash_id":value.trash_id,"operation_id":value.operation_id,"source":value.canonical_source,"consumer_generation":value.consumer,"purge_after_ms":value.purge_after_ms}))?;
            }
            Ok(())
        })
    }
    pub(super) fn mark_file_purge(&self, rows: &[Stored]) -> Result<()> {
        for row in rows {
            self.verify_trash_directory(&row.descriptor.trash_id)?;
            self.verify_trash_file(
                &self.file_trash_path(&row.descriptor.trash_id, false)?,
                &row.descriptor,
            )?;
        }
        self.atomic(|db|db.write(|conn| {for row in rows {conn.execute("UPDATE file_trash SET state='PURGING' WHERE trash_id=?1 AND state='TRASHED'",[&row.descriptor.trash_id])?;}Ok(())}))
    }
    pub(super) fn recover_file_trash(&self) -> Result<()> {
        let rows = self.file_trash_rows()?;
        let mut restoring = BTreeMap::<String, Vec<Descriptor>>::new();
        for row in rows {
            let value = &row.descriptor;
            match row.state.as_str() {
                "MOVING" => {
                    self.move_cleanup_file(value, false)?;
                    self.write(|conn|{conn.execute("UPDATE file_trash SET state='TRASHED' WHERE trash_id=?1 AND state='MOVING'",[&value.trash_id])?;Ok(())})?;
                }
                "RESTORING" => {
                    restoring
                        .entry(value.operation_id.clone())
                        .or_default()
                        .push(value.clone());
                }
                "PURGING" => {
                    self.verify_trash_directory(&value.trash_id)?;
                    let data = self.file_trash_path(&value.trash_id, false)?;
                    if data.exists() {
                        self.verify_trash_file(&data, value)?;
                        fs::remove_file(&data)?;
                        File::open(data.parent().unwrap())?.sync_all()?;
                    }
                    self.write(|conn|{conn.execute("UPDATE file_trash SET state='PURGED' WHERE trash_id=?1 AND state='PURGING'",[&value.trash_id])?;Ok(())})?;
                }
                "TRASHED" | "RESTORED" | "PURGED" => {}
                _ => bail!("unsupported file trash transition"),
            }
        }
        for (_, values) in restoring {
            for value in &values {
                self.move_cleanup_file(value, true)?;
            }
            let register = self.restore_registrations_allowed(&values)?;
            self.atomic(|db| {
                if register {for value in &values {db.record_artifact(&value.record)?;}}
                for value in &values {
                    db.write(|conn|{conn.execute("UPDATE file_trash SET state='RESTORED' WHERE trash_id=?1 AND state='RESTORING'",[&value.trash_id])?;Ok(())})?;
                    db.append_event("cleanup",Some(&value.run),"input file restored",serde_json::json!({"trash_id":value.trash_id,"operation_id":value.operation_id,"registered":register}))?;
                }
                Ok(())
            })?;
        }
        Ok(())
    }
}
