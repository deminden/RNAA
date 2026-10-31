//! Revalidate the recorded QC/count admission and selected reference identities.
use crate::quality;
use anyhow::{Context, Result, bail};
use rnaa_core::{
    Database, ProjectConfig,
    generations::GenerationHead,
    references::{ReferenceScope, ReferenceSelection},
};

pub(super) struct Guard {
    eligibility: quality::Snapshot,
    reference: ReferenceSelection,
    eligibility_digest: String,
    reference_digest: String,
}

impl Guard {
    pub(super) fn eligible(&self) -> &std::collections::BTreeSet<String> {
        &self.eligibility.eligible
    }
    pub(super) fn reference(&self) -> &ReferenceSelection {
        &self.reference
    }
    pub(super) fn bind<P>(
        &self,
        request: &mut rnaa_stages::analysis_request::AnalysisRequest<P>,
    ) -> Result<()> {
        self.eligibility.bind(request)?;
        crate::reference::bind(request, &self.reference)
    }
    pub(super) fn capture(
        db: &Database,
        config: &ProjectConfig,
        heads: &[GenerationHead],
    ) -> Result<Self> {
        let eligibility = quality::Snapshot::observe(db, config)?;
        let eligibility_digest = eligibility.digest()?;
        let reference = db
            .selected_reference()?
            .context("no selected reference; cannot check current admission")?;
        let reference_digest = reference.fingerprint(ReferenceScope::Analysis)?;
        for head in heads {
            let request = db.scientific_request(&head.request_sha256)?.context("selected generation lacks a verified scientific request; cannot check current admission")?;
            if request.dependencies.get("eligibility") != Some(&eligibility_digest) {
                bail!(
                    "selected {} generation does not match current QC eligibility or quantified abundance inputs",
                    head.output_group.as_str()
                );
            }
            if request.dependencies.get("reference") != Some(&reference_digest) {
                bail!(
                    "selected {} generation does not match the current selected reference",
                    head.output_group.as_str()
                );
            }
        }
        let guard = Self {
            eligibility,
            reference,
            eligibility_digest,
            reference_digest,
        };
        guard.verify_current(db, config)?;
        Ok(guard)
    }

    pub(super) fn verify_current(&self, db: &Database, config: &ProjectConfig) -> Result<()> {
        self.eligibility.verify_current(db, config)?;
        self.reference
            .verify_current(db, config, ReferenceScope::Analysis)?;
        Ok(())
    }

    pub(super) fn verify_catalog(&self, db: &Database, config: &ProjectConfig) -> Result<()> {
        self.eligibility.verify_catalog(db, config)?;
        self.reference
            .verify_catalog(db, config, ReferenceScope::Analysis)?;
        Ok(())
    }

    pub(super) fn proof(&self) -> serde_json::Value {
        serde_json::json!({"eligibility_sha256":self.eligibility_digest,"reference_sha256":self.reference_digest,"full_freshness":false})
    }
}
