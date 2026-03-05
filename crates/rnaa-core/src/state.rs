use std::fmt::{Display, Formatter};
use std::str::FromStr;

use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RunState {
    New,
    Resolved,
    Downloading,
    Downloaded,
    Verified,
    DownloadFailed,
    VerifyFailed,
    QuantRunning,
    QuantDone,
    QuantFailed,
    DeRunning,
    DeDone,
    DeFailed,
    CorrRunning,
    CorrDone,
    CorrFailed,
    Cleaned,
}

impl RunState {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::New => "NEW",
            Self::Resolved => "RESOLVED",
            Self::Downloading => "DOWNLOADING",
            Self::Downloaded => "DOWNLOADED",
            Self::Verified => "VERIFIED",
            Self::DownloadFailed => "DOWNLOAD_FAILED",
            Self::VerifyFailed => "VERIFY_FAILED",
            Self::QuantRunning => "QUANT_RUNNING",
            Self::QuantDone => "QUANT_DONE",
            Self::QuantFailed => "QUANT_FAILED",
            Self::DeRunning => "DE_RUNNING",
            Self::DeDone => "DE_DONE",
            Self::DeFailed => "DE_FAILED",
            Self::CorrRunning => "CORR_RUNNING",
            Self::CorrDone => "CORR_DONE",
            Self::CorrFailed => "CORR_FAILED",
            Self::Cleaned => "CLEANED",
        }
    }
}

impl Display for RunState {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for RunState {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "NEW" => Ok(Self::New),
            "RESOLVED" => Ok(Self::Resolved),
            "DOWNLOADING" => Ok(Self::Downloading),
            "DOWNLOADED" => Ok(Self::Downloaded),
            "VERIFIED" => Ok(Self::Verified),
            "DOWNLOAD_FAILED" => Ok(Self::DownloadFailed),
            "VERIFY_FAILED" => Ok(Self::VerifyFailed),
            "QUANT_RUNNING" => Ok(Self::QuantRunning),
            "QUANT_DONE" => Ok(Self::QuantDone),
            "QUANT_FAILED" => Ok(Self::QuantFailed),
            "DE_RUNNING" => Ok(Self::DeRunning),
            "DE_DONE" => Ok(Self::DeDone),
            "DE_FAILED" => Ok(Self::DeFailed),
            "CORR_RUNNING" => Ok(Self::CorrRunning),
            "CORR_DONE" => Ok(Self::CorrDone),
            "CORR_FAILED" => Ok(Self::CorrFailed),
            "CLEANED" => Ok(Self::Cleaned),
            _ => bail!("unknown run state: {value}"),
        }
    }
}
