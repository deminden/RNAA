use super::*;
#[test]
fn failed_upgrade_rolls_back_schema_and_version_records() -> Result<()> {
    let conn = Connection::open_in_memory()?;
    conn.execute_batch("CREATE TABLE artifacts (project_id TEXT, run_accession TEXT, kind TEXT, path TEXT, blob_id TEXT)")?;
    assert!(Database::migrate(&conn).is_err());
    let count: i64 = conn.query_row("SELECT count(*) FROM sqlite_master WHERE name IN ('runs', 'schema_migrations', 'retained_inputs')", [], |row| row.get(0))?;
    assert_eq!(count, 0);
    Ok(())
}
