//! Catalog ownership for export directories. Filesystem reconciliation holds the
//! project analysis lease; these records alone never authorize file deletion.
use crate::Database;
use anyhow::{Context, Result, bail};
use rusqlite::params;
use std::path::{Component, Path};

#[derive(Debug)]
pub struct PendingExport {
    pub attempt_id: String,
    pub stage: String,
    pub staged_name: String,
    pub device: u64,
    pub inode: u64,
}

fn valid_name(stage: &str, attempt: &str, name: &str) -> bool {
    matches!(stage, "export" | "survey")
        && uuid::Uuid::parse_str(attempt).is_ok()
        && name.starts_with(&format!(".rnaa-{stage}-{attempt}-"))
        && name.len() > format!(".rnaa-{stage}-{attempt}-").len()
        && Path::new(name).file_name().and_then(|part| part.to_str()) == Some(name)
        && matches!(
            Path::new(name).components().collect::<Vec<_>>().as_slice(),
            [Component::Normal(_)]
        )
}

impl Database {
    fn check_export_lease(&self, conn: &rusqlite::Connection) -> Result<()> {
        let identity = self
            .job
            .as_ref()
            .context("export recovery requires an owned job")?;
        identity.check(conn)?;
        let valid: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM job_attempts WHERE attempt_id=?1 AND project_id=?2 AND resource='analysis' AND stage IN ('export','survey'))",params![identity.attempt_id,self.project_id()],|row|row.get(0))?;
        if !valid {
            bail!("export recovery requires an export or survey analysis lease");
        }
        Ok(())
    }

    pub fn register_export_stage(
        &self,
        stage: &str,
        name: &str,
        device: u64,
        inode: u64,
    ) -> Result<()> {
        let identity = self
            .job
            .as_ref()
            .context("export registration requires an owned job")?;
        if !valid_name(stage, &identity.attempt_id, name) {
            bail!("invalid export staging name");
        }
        let device =
            i64::try_from(device).context("export device identity exceeds catalog range")?;
        let inode = i64::try_from(inode).context("export inode identity exceeds catalog range")?;
        self.write(|conn| {
            self.check_export_lease(conn)?;
            let actual: String = conn.query_row("SELECT stage FROM job_attempts WHERE attempt_id=?1 AND project_id=?2", params![identity.attempt_id, self.project_id()], |row| row.get(0))?;
            if actual != stage { bail!("export stage does not match owned job"); }
            conn.execute("INSERT INTO export_publications(attempt_id, project_id, stage, staged_name, device, inode, state) VALUES (?1,?2,?3,?4,?5,?6,'STAGED')", params![identity.attempt_id,self.project_id(),stage,name,device,inode])?;
            Ok(())
        })
    }

    pub fn mark_export_published(&self) -> Result<()> {
        let identity = self
            .job
            .as_ref()
            .context("export publication requires an owned job")?;
        self.write(|conn| {
            self.check_export_lease(conn)?;
            if conn.execute("UPDATE export_publications SET state='PUBLISHED' WHERE attempt_id=?1 AND project_id=?2 AND state='STAGED'",params![identity.attempt_id,self.project_id()])? != 1 { bail!("export publication record is missing or changed"); }
            Ok(())
        })
    }

    pub fn pending_export_recovery(&self) -> Result<Vec<PendingExport>> {
        if self.job.is_none() {
            bail!("export recovery requires the project analysis lease");
        }
        self.read(|conn| {
            self.check_export_lease(conn)?;
            let mut query=conn.prepare("SELECT p.attempt_id,p.stage,p.staged_name,p.device,p.inode FROM export_publications p JOIN job_attempts a ON a.attempt_id=p.attempt_id AND a.project_id=p.project_id AND a.stage=p.stage WHERE p.project_id=?1 AND p.state='STAGED' AND a.state<>'RUNNING' ORDER BY p.attempt_id")?;
            let rows=query.query_map([self.project_id()], |row| Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?,row.get::<_,i64>(3)?,row.get::<_,i64>(4)?)))?;
            let mut records=Vec::new();
            for row in rows { let (attempt_id,stage,staged_name,device,inode)=row?;
                if !valid_name(&stage,&attempt_id,&staged_name) { bail!("invalid export recovery record"); }
                records.push(PendingExport { attempt_id,stage,staged_name,device:u64::try_from(device)?,inode:u64::try_from(inode)? });
            }
            Ok(records)
        })
    }

    pub fn finish_export_recovery(&self, attempt: &str, state: &str) -> Result<()> {
        if self.job.is_none() || !matches!(state, "RECOVERED" | "ORPHANED" | "MISSING") {
            bail!("invalid export recovery transition");
        }
        self.write(|conn| {
            self.check_export_lease(conn)?;
            let (stage, name): (String, String) = conn.query_row("SELECT stage, staged_name FROM export_publications WHERE project_id=?1 AND attempt_id=?2",params![self.project_id(),attempt],|row| Ok((row.get(0)?,row.get(1)?)))?;
            if !valid_name(&stage,attempt,&name) { bail!("invalid export recovery record"); }
            if conn.execute("UPDATE export_publications SET state=?1 WHERE project_id=?2 AND attempt_id=?3 AND state='STAGED' AND EXISTS (SELECT 1 FROM job_attempts a WHERE a.attempt_id=export_publications.attempt_id AND a.project_id=export_publications.project_id AND a.stage=export_publications.stage AND a.state<>'RUNNING')",params![state,self.project_id(),attempt])? != 1 { bail!("export recovery record changed"); }
            let path = match state {
                "RECOVERED" => Some(format!(".rnaa-recovered-{}", &name[1..])),
                "ORPHANED" => Some(name[1..].to_owned()),
                _ => None,
            };
            conn.execute("INSERT INTO events(project_id,ts,stage,run_accession,message,context_json) VALUES (?1,?2,?3,NULL,'interrupted export reconciled',?4)",params![self.project_id(),crate::util::now_rfc3339(),stage,serde_json::json!({"attempt_id":attempt,"state":state,"exports_relative_path":path}).to_string()])?;
            Ok(())
        })
    }
}
