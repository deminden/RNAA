//! Completion predicates for legacy run progress; QC eligibility remains independent.
use rnaa_core::state::RunState;

pub(crate) fn is_quant_done_state(state: RunState) -> bool {
    matches!(
        state,
        RunState::QuantDone
            | RunState::DeRunning
            | RunState::DeDone
            | RunState::DeFailed
            | RunState::CorrRunning
            | RunState::CorrDone
            | RunState::CorrFailed
            | RunState::Cleaned
    )
}

pub(crate) fn is_preprocess_done_state(state: RunState) -> bool {
    matches!(
        state,
        RunState::PreprocessDone
            | RunState::QuantRunning
            | RunState::QuantDone
            | RunState::DeRunning
            | RunState::DeDone
            | RunState::DeFailed
            | RunState::CorrRunning
            | RunState::CorrDone
            | RunState::CorrFailed
            | RunState::Cleaned
    )
}

pub(crate) fn is_download_complete_state(state: RunState) -> bool {
    matches!(
        state,
        RunState::Downloaded
            | RunState::Verified
            | RunState::PreprocessRunning
            | RunState::PreprocessDone
            | RunState::PreprocessFailed
            | RunState::QuantRunning
            | RunState::QuantDone
            | RunState::QuantFailed
            | RunState::DeRunning
            | RunState::DeDone
            | RunState::DeFailed
            | RunState::CorrRunning
            | RunState::CorrDone
            | RunState::CorrFailed
            | RunState::Cleaned
    )
}
