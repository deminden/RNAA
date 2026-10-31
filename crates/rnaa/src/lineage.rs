//! Read-only, machine-readable provenance inspection with optional payload checks.
use crate::cli::LineageArgs;
use anyhow::{Context, Result, bail};
use rnaa_core::Database;
use std::collections::{BTreeMap, VecDeque};

pub(crate) fn inspect(db: &Database, args: LineageArgs) -> Result<()> {
    let _pin =
        rnaa_core::generations::GenerationStorePin::acquire(db.path(), "lineage inspection")?;
    let heads = db.generation_heads()?;
    let selected_ids = heads
        .iter()
        .map(|head| head.generation_id.clone())
        .collect::<std::collections::BTreeSet<_>>();
    let mut queue = heads
        .into_iter()
        .filter(|head| {
            args.group
                .as_ref()
                .is_none_or(|group| *group == head.output_group.as_str())
                && args
                    .subject
                    .as_ref()
                    .is_none_or(|subject| *subject == head.subject)
        })
        .collect::<VecDeque<_>>();
    if args.all_requests {
        let selected = queue.iter().cloned().collect::<Vec<_>>();
        for head in selected {
            queue.extend(db.committed_generation_history(head.output_group, &head.subject)?);
        }
    }
    let mut entries = BTreeMap::new();
    let mut failed = false;
    while let Some(head) = queue.pop_front() {
        if entries.contains_key(&head.generation_id) {
            continue;
        }
        let lineage = db.generation_lineage(&head)?;
        let input_binding = if selected_ids.contains(&head.generation_id) {
            db.generation_input_binding(head.output_group, &head.subject)?
        } else {
            None
        };
        if args.history {
            for parent in lineage.dependencies.iter().chain(
                input_binding
                    .iter()
                    .flat_map(|binding| &binding.dependencies),
            ) {
                queue.push_back(
                    db.committed_generation_head(&parent.producer_generation_id)?
                        .context(
                            "recorded producer is missing from this project's committed history",
                        )?,
                );
            }
        }
        let verified = args
            .verify
            .then(|| db.verify_generation_head(&head))
            .transpose()?;
        failed |= verified == Some(false);
        let binding_verified = input_binding.as_ref().and_then(|binding| {
            args.verify
                .then(|| db.verify_input_binding(binding).is_ok())
        });
        failed |= binding_verified == Some(false);
        entries.insert(
            head.generation_id.clone(),
            serde_json::json!({"lineage":lineage,"payload_verified":verified,"selected":selected_ids.contains(&head.generation_id),"input_binding":input_binding,"input_binding_verified":binding_verified}),
        );
    }
    println!(
        "{}",
        serde_json::to_string_pretty(
            &serde_json::json!({"schema_version":1,"generations":entries.into_values().collect::<Vec<_>>()})
        )?
    );
    if failed {
        bail!("one or more generation payloads failed verification; inspect the JSON report");
    }
    Ok(())
}
