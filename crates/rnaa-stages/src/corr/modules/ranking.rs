//! Explicit finite t-score policy; undefined genes are omitted from the universe.
pub(super) fn t_score(r: f64, samples: usize) -> f64 {
    let r = r.clamp(-0.999_999_9, 0.999_999_9);
    r * (((samples - 2) as f64) / (1.0 - r * r)).sqrt()
}
