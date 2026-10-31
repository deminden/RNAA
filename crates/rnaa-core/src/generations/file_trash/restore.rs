use super::*;

impl Database {
    pub(super) fn restore_registrations_allowed(&self, values: &[Descriptor]) -> Result<bool> {
        for value in values {
            let run = self
                .get_run(&value.run)?
                .context("file restore run disappeared")?;
            if digest(serde_json::to_vec(&run.remote_files)?) != value.source_selection
                || run.library_layout != value.layout
                || canonical_missing_path(Path::new(&value.record.path))?
                    != Path::new(&value.canonical_source)
            {
                return Ok(false);
            }
            if let Some(id) = &value.producer {
                let producer = self
                    .committed_seal(id)?
                    .context("file restore producer disappeared")?;
                if self
                    .generation_head(producer.head.output_group, &value.run)?
                    .is_none_or(|head| head.generation_id != *id)
                {
                    return Ok(false);
                }
            }
            let family = |kind: ArtifactKind| match kind {
                ArtifactKind::FastqR1 | ArtifactKind::FastqR2 | ArtifactKind::FastqSingle => 0,
                ArtifactKind::FastqTrimmedR1
                | ArtifactKind::FastqTrimmedR2
                | ArtifactKind::FastqTrimmedSingle => 1,
                ArtifactKind::Sra => 2,
                _ => 3,
            };
            for active in self.list_artifacts_for_run(Some(&value.run))? {
                if family(active.kind) == family(value.record.kind)
                    && !values
                        .iter()
                        .any(|expected| expected.record.path == active.path)
                {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
    pub fn restore_file_trash(&self, id: &str) -> Result<serde_json::Value> {
        let _lock = super::super::pins::maintenance_lock(self)?
            .context("project inputs are in use; retry restore")?;
        let mut db = self.clone();
        db.generation_maintenance = true;
        db.recover_file_trash()?;
        db.restore_file_operation(id)
    }
    fn restore_file_operation(&self, id: &str) -> Result<serde_json::Value> {
        let rows = self.file_trash_rows()?;
        let selected = rows
            .iter()
            .find(|row| row.descriptor.trash_id == id)
            .context("unknown file trash identity")?;
        let operation = &selected.descriptor.operation_id;
        let group = rows
            .iter()
            .filter(|row| row.descriptor.operation_id == *operation)
            .cloned()
            .collect::<Vec<_>>();
        if group.iter().all(|row| row.state == "RESTORED") {
            return Ok(
                serde_json::json!({"operation_id":operation,"restored":group.iter().map(|row|&row.descriptor.trash_id).collect::<Vec<_>>(),"already_restored":true}),
            );
        }
        if selected.state != "TRASHED" {
            bail!("only retained file trash can be restored");
        }
        if group.iter().any(|row| row.state != "TRASHED") {
            bail!("the complete cleanup operation is no longer restorable");
        }
        for row in &group {
            if self.copy_staging_path(&row.descriptor, true)?.exists() {
                bail!("restore copy staging identity already exists");
            }
            self.verify_trash_directory(&row.descriptor.trash_id)?;
            self.verify_trash_file(
                &self.file_trash_path(&row.descriptor.trash_id, false)?,
                &row.descriptor,
            )?;
        }
        self.atomic(|db|db.write(|conn| {for row in &group {conn.execute("UPDATE file_trash SET state='RESTORING' WHERE trash_id=?1 AND state='TRASHED'",[&row.descriptor.trash_id])?;}Ok(())}))?;
        self.recover_file_trash()?;
        Ok(
            serde_json::json!({"operation_id":operation,"restored":group.iter().map(|row|&row.descriptor.trash_id).collect::<Vec<_>>()}),
        )
    }
    fn restorable_file_groups(&self, run: &str) -> Result<Vec<Vec<Descriptor>>> {
        let rows = self.file_trash_rows_filtered(Some(run), None)?;
        let mut groups = BTreeMap::<String, Vec<&Stored>>::new();
        for row in &rows {
            groups
                .entry(row.descriptor.operation_id.clone())
                .or_default()
                .push(row);
        }
        let mut result = Vec::new();
        for group in groups.values() {
            if !group.iter().any(|row| fastq(row.descriptor.record.kind))
                || group
                    .iter()
                    .any(|row| !matches!(row.state.as_str(), "MOVING" | "TRASHED"))
            {
                continue;
            }
            let values = group
                .iter()
                .map(|row| row.descriptor.clone())
                .collect::<Vec<_>>();
            if !self.restore_registrations_allowed(&values)? {
                continue;
            }
            let mut valid = true;
            for value in &values {
                self.verify_trash_directory(&value.trash_id)?;
                let trash = self.file_trash_path(&value.trash_id, false)?;
                let original = self.original_file_path(value, false)?;
                if (trash.exists() && self.verify_trash_file(&trash, value).is_err())
                    || (original.exists() && self.verify_trash_file(&original, value).is_err())
                    || (!trash.exists() && !original.exists())
                {
                    valid = false;
                    break;
                }
            }
            if valid {
                result.push(values);
            }
        }
        Ok(result)
    }
    pub fn restorable_file_trash(&self, run: &str) -> Result<bool> {
        let _pin = if self.generation_maintenance {
            None
        } else {
            Some(super::super::pins::read_lock(self.path())?)
        };
        Ok(!self.restorable_file_groups(run)?.is_empty())
    }
    pub fn restore_run_file_trash(&self, run: &str) -> Result<bool> {
        let _lock = super::super::pins::maintenance_lock(self)?
            .context("project inputs are in use; retry recovery")?;
        let mut db = self.clone();
        db.generation_maintenance = true;
        db.recover_file_trash()?;
        if let Some(values) = db.restorable_file_groups(run)?.first() {
            db.restore_file_operation(&values[0].trash_id)?;
            return Ok(true);
        }
        Ok(false)
    }
}
