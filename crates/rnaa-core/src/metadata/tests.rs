use super::*;

#[test]
fn numeric_looking_conditions_are_factors_with_declared_reference_order() -> Result<()> {
    let mut schema = MetadataConfig::default();
    schema.columns.insert(
        "condition".into(),
        ColumnSpec {
            kind: ColumnType::Factor,
            levels: vec!["2".into(), "10".into(), "1".into()],
            reference: Some("2".into()),
        },
    );
    let headers = vec!["condition".into()];
    let rows = vec![
        vec!["1".into()],
        vec!["2".into()],
        vec!["1".into()],
        vec!["2".into()],
    ];
    let frozen = FrozenMetadata::resolve(&headers, &rows, "~ condition", &schema)?;
    let FrozenColumn::Factor {
        levels, reference, ..
    } = &frozen.columns["condition"]
    else {
        panic!("numeric-looking condition must remain categorical")
    };
    assert_eq!(levels, &["2", "1"]);
    assert_eq!(reference, "2");
    assert_eq!(
        frozen,
        serde_json::from_slice::<FrozenMetadata>(&serde_json::to_vec(&frozen)?)?
    );
    assert!(
        FrozenMetadata::resolve(&headers, &[vec!["1".into()]], "~ condition", &schema).is_err()
    );
    Ok(())
}

#[test]
fn explicit_numeric_roles_missing_rules_and_forbidden_identifiers_are_checked() -> Result<()> {
    let mut schema = MetadataConfig::default();
    assert!(schema.normalize("age", "21").is_err());
    schema
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    assert_eq!(schema.normalize("age", " 21.0 ")?, "21");
    for value in ["NaN", "inf", "young", ""] {
        assert!(schema.normalize("age", value).is_err());
    }
    schema.missing_values.push("NA".into());
    assert!(schema.normalize("condition", "NA").is_err());
    schema.columns.insert(
        "source file (URL)".into(),
        ColumnSpec::new(ColumnType::Ignore),
    );
    schema.validate()?;
    schema.columns.insert(
        "sample_accession".into(),
        ColumnSpec::new(ColumnType::Factor),
    );
    assert!(schema.validate().is_err());
    Ok(())
}

#[test]
fn formulas_discover_real_covariates_without_confusing_transform_arguments() -> Result<()> {
    let mut schema = MetadataConfig::default();
    schema
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    let headers = vec![
        "age".into(),
        "condition".into(),
        "run_accession".into(),
        "unused".into(),
    ];
    let expected = BTreeSet::from(["age".to_string(), "condition".to_string()]);
    assert_eq!(
        variables(
            "~ relevel(condition, ref='B') * I(age^2) + poly(age, degree=2, raw=TRUE)",
            &headers,
            &schema
        )?,
        expected
    );
    assert_eq!(variables("~ .", &headers, &schema)?, expected);
    assert_eq!(
        variables(
            "~ relevel(condition, ref='B~C') + I(age * .5e-3)",
            &headers,
            &schema
        )?,
        expected
    );
    assert!(variables("~ unknown", &headers, &schema).is_err());
    assert!(variables("~ run_accession", &headers, &schema).is_err());
    assert!(variables("~ relevel(condition, ref='B)", &headers, &schema).is_err());
    Ok(())
}

#[test]
fn malformed_schema_fields_cannot_be_silently_ignored() -> Result<()> {
    assert!(
        toml::from_str::<MetadataConfig>("minimum_replicates = 0")?
            .validate()
            .is_err()
    );
    assert!(toml::from_str::<MetadataConfig>("min_replicates = 2").is_err());
    assert!(
        toml::from_str::<ColumnSpec>("type='numeric'\nreference='1'")?.kind == ColumnType::Numeric
    );
    let mut schema = MetadataConfig::default();
    schema.columns.insert(
        "age".into(),
        toml::from_str("type='numeric'\nreference='1'")?,
    );
    assert!(schema.validate().is_err());
    assert!(toml::from_str::<ColumnSpec>("type='factor'\nreferance='A'").is_err());
    Ok(())
}
