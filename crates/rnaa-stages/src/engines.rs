//! Resolved engine identities, generated from the workspace lockfile at build time.
pub const KALLISTORS: &str = env!("RNAA_KALLISTORS_VERSION");
pub const RSDESEQ2: &str = env!("RNAA_RSDESEQ2_VERSION");
pub const RSFGSEA: &str = env!("RNAA_RSFGSEA_VERSION");
pub const MINCORR: &str = env!("RNAA_MINCORR_VERSION");
pub const MINCORR_REVISION: &str = env!("RNAA_MINCORR_REVISION");
pub const DEPENDENCY_LOCK_SHA256: &str = env!("RNAA_DEPENDENCY_LOCK_SHA256");
