//! A permit follows each batch into its ordered result, covering every queue.
use anyhow::{bail, Result};
use std::sync::{Arc, Condvar, Mutex};

#[derive(Default)]
struct State {
    active: usize,
    peak: usize,
    cancelled: bool,
}
pub(super) struct BatchWindow {
    limit: usize,
    state: Mutex<State>,
    changed: Condvar,
}
pub(crate) struct BatchPermit {
    window: Arc<BatchWindow>,
}
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub(super) struct WindowStats {
    pub limit: usize,
    pub peak: usize,
    pub active: usize,
}
pub(super) struct PanicGuard(pub Arc<BatchWindow>);
pub(super) struct ScopeGuard(pub Arc<BatchWindow>);
impl Drop for ScopeGuard {
    fn drop(&mut self) {
        self.0.cancel();
    }
}
impl Drop for PanicGuard {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.cancel();
        }
    }
}
impl BatchWindow {
    pub fn new(limit: usize) -> Arc<Self> {
        Arc::new(Self {
            limit,
            state: Mutex::new(State::default()),
            changed: Condvar::new(),
        })
    }
    pub fn acquire(self: &Arc<Self>) -> Result<Arc<BatchPermit>> {
        let mut state = self
            .state
            .lock()
            .map_err(|_| anyhow::anyhow!("native batch window is poisoned"))?;
        while state.active >= self.limit && !state.cancelled {
            state = self
                .changed
                .wait(state)
                .map_err(|_| anyhow::anyhow!("native batch window is poisoned"))?;
        }
        if state.cancelled {
            bail!("native preprocessing pipeline was cancelled");
        }
        state.active += 1;
        state.peak = state.peak.max(state.active);
        Ok(Arc::new(BatchPermit {
            window: Arc::clone(self),
        }))
    }
    pub fn cancel(&self) {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        state.cancelled = true;
        self.changed.notify_all();
    }
    pub fn stats(&self) -> Result<WindowStats> {
        let state = self
            .state
            .lock()
            .map_err(|_| anyhow::anyhow!("native batch window is poisoned"))?;
        Ok(WindowStats {
            limit: self.limit,
            peak: state.peak,
            active: state.active,
        })
    }
}
impl Drop for BatchPermit {
    fn drop(&mut self) {
        let mut state = self
            .window
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if state.active > 0 {
            state.active -= 1;
        } else {
            state.cancelled = true;
        }
        self.window.changed.notify_all();
    }
}
