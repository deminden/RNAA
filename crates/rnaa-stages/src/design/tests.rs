use super::*;
use rnaa_core::metadata::{ColumnSpec, ColumnType};

fn inputs() -> (Vec<String>, Vec<Vec<String>>, Vec<String>) {
    (
        vec!["condition".into(), "age".into()],
        (0..8)
            .map(|i| {
                vec![
                    if i < 4 { "1".into() } else { "2".into() },
                    (20 + i).to_string(),
                ]
            })
            .collect(),
        (0..8).map(|i| format!("S{i}")).collect(),
    )
}
#[test]
fn explicit_numeric_and_factor_contracts_share_the_exact_model_and_reference() -> Result<()> {
    let (headers, rows, ids) = inputs();
    let mut schema = MetadataConfig::default();
    schema
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    schema.columns.insert(
        "condition".into(),
        ColumnSpec {
            kind: ColumnType::Factor,
            levels: vec!["2".into(), "1".into()],
            reference: Some("2".into()),
        },
    );
    let contrast = ContrastSpec {
        name: "effect".into(),
        factor: "condition".into(),
        level_a: "1".into(),
        level_b: "2".into(),
    };
    let de = PreparedDesign::new(
        &headers,
        &rows,
        &ids,
        "~ age + condition",
        &schema,
        &[contrast],
        Purpose::DifferentialExpression,
    )?;
    let corr = PreparedDesign::new(
        &headers,
        &rows,
        &ids,
        "~ age + condition",
        &schema,
        &[],
        Purpose::Correlation,
    )?;
    assert_eq!(de.contract.matrix, corr.contract.matrix);
    assert_eq!(de.contract.factor_references["condition"], "2");
    assert_eq!(de.contract.rank, 3);
    assert_eq!(de.contract.residual_degrees_of_freedom, 5);
    assert_eq!(de.contract.contrasts[0].weights.len(), 3);
    let intercept = PreparedDesign::new(
        &headers,
        &rows,
        &ids,
        "~ 1",
        &schema,
        &[],
        Purpose::DifferentialExpression,
    )?;
    assert_eq!(intercept.contract.matrix, vec![1.0; 8]);
    Ok(())
}
#[test]
fn rank_replication_and_contrast_fail_before_expression_fitting() -> Result<()> {
    let (headers, mut rows, ids) = inputs();
    let mut schema = MetadataConfig::default();
    schema
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    for row in &mut rows {
        row[1] = row[0].clone();
    }
    assert!(
        PreparedDesign::new(
            &headers,
            &rows,
            &ids,
            "~ age + condition",
            &schema,
            &[],
            Purpose::DifferentialExpression
        )
        .unwrap_err()
        .to_string()
        .contains("rank deficient")
    );
    let contrast = ContrastSpec {
        name: "effect".into(),
        factor: "condition".into(),
        level_a: "1".into(),
        level_b: "2".into(),
    };
    let few = vec![
        vec!["1".into(), "20".into()],
        vec!["2".into(), "21".into()],
        vec!["2".into(), "22".into()],
    ];
    assert!(
        PreparedDesign::new(
            &headers,
            &few,
            &ids[..3],
            "~ condition",
            &schema,
            std::slice::from_ref(&contrast),
            Purpose::DifferentialExpression
        )
        .unwrap_err()
        .to_string()
        .contains("minimum_replicates")
    );
    let mut missing = contrast;
    missing.level_a = "absent".into();
    assert!(
        PreparedDesign::new(
            &headers,
            &rows,
            &ids,
            "~ condition",
            &schema,
            &[missing],
            Purpose::DifferentialExpression
        )
        .is_err()
    );
    assert!(
        PreparedDesign::new(
            &headers,
            &few[..2],
            &ids[..2],
            "~ condition",
            &schema,
            &[],
            Purpose::DifferentialExpression
        )
        .unwrap_err()
        .to_string()
        .contains("saturated")
    );
    Ok(())
}
#[test]
fn transforms_and_dot_use_declared_roles_and_offsets_are_not_silently_ignored() -> Result<()> {
    let (headers, rows, ids) = inputs();
    let mut schema = MetadataConfig::default();
    schema
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    let squared = PreparedDesign::new(
        &headers,
        &rows,
        &ids,
        "~ I(age^2) + condition",
        &schema,
        &[],
        Purpose::Correlation,
    )?;
    assert_eq!(squared.contract.rank, 3);
    assert!(
        PreparedDesign::new(
            &headers,
            &rows,
            &ids,
            "~ condition + offset(age)",
            &schema,
            &[],
            Purpose::Correlation
        )
        .is_err()
    );
    let offset = PreparedDesign::new(
        &headers,
        &rows,
        &ids,
        "~ condition + offset(age)",
        &schema,
        &[],
        Purpose::DifferentialExpression,
    )?;
    assert_eq!(
        offset.contract.offsets,
        (20..28).map(f64::from).collect::<Vec<_>>()
    );
    assert_eq!(
        PreparedDesign::new(
            &headers,
            &rows,
            &ids,
            "~ .",
            &schema,
            &[],
            Purpose::Correlation
        )?
        .contract
        .matrix
        .len(),
        24
    );
    let contrast = ContrastSpec {
        name: "effect".into(),
        factor: "condition".into(),
        level_a: "1".into(),
        level_b: "2".into(),
    };
    let releveled = PreparedDesign::new(
        &headers,
        &rows,
        &ids,
        "~ relevel(condition, ref='2')",
        &schema,
        &[contrast],
        Purpose::DifferentialExpression,
    )?;
    assert_eq!(
        releveled.contract.contrasts[0].reference.as_deref(),
        Some("2")
    );
    assert_eq!(releveled.contract.contrasts[0].weights, vec![0.0, 1.0]);
    let contrast = ContrastSpec {
        name: "effect".into(),
        factor: "condition".into(),
        level_a: "1".into(),
        level_b: "2".into(),
    };
    let indicators = PreparedDesign::new(
        &headers,
        &rows,
        &ids,
        "~ 0 + condition",
        &schema,
        &[contrast],
        Purpose::DifferentialExpression,
    )?;
    assert_eq!(indicators.contract.contrasts[0].reference, None);
    assert_eq!(indicators.contract.contrasts[0].weights, vec![1.0, -1.0]);
    Ok(())
}
