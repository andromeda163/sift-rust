//! Fast exponential approximation (exact port of VLFeat's fast_expn from sift.c).

use std::sync::OnceLock;

const EXPN_SZ: usize = 256;
const EXPN_MAX: f64 = 25.0;

pub(crate) static EXPN_TAB: OnceLock<[f64; EXPN_SZ + 1]> = OnceLock::new();

fn get_expn_tab() -> &'static [f64; EXPN_SZ + 1] {
    EXPN_TAB.get_or_init(|| {
        let mut tab = [0.0f64; EXPN_SZ + 1];
        for k in 0..=EXPN_SZ {
            tab[k] = (-(k as f64) * (EXPN_MAX / EXPN_SZ as f64)).exp();
        }
        tab
    })
}

/// Initialize the fast_expn lookup table (safe to call multiple times; thread-safe).
pub fn fast_expn_init() {
    let _ = get_expn_tab();
}

/// Fast exp(-x) approximation using table + linear interpolation.
/// Thread-safe: uses OnceLock (no static mut/unsafe).
/// OOB-safe: `if x >= EXPN_MAX` returns 0.0 before any indexing;
/// index guard `i < EXPN_SZ` before `[i+1]` as defense-in-depth.
/// Table ends at exp(-25.0); boundary returns 0.0 (safe, no UB/panic).
pub fn fast_expn(x: f64) -> f64 {
    if x >= EXPN_MAX {
        return 0.0;
    }
    let tab = get_expn_tab();
    let xs = x * (EXPN_SZ as f64 / EXPN_MAX);
    let i = xs.floor() as usize;
    let r = xs - i as f64;
    let a = tab[i];
    // Guard: i+1 is valid only if i < EXPN_SZ (257 elements, indices 0..256)
    let b = if i < EXPN_SZ { tab[i + 1] } else { a };
    a + r * (b - a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_fast_expn_init_and_values() {
        fast_expn_init();
        let tab = get_expn_tab();
        assert!((tab[0] - 1.0).abs() < 1e-10);
        assert!((fast_expn(0.0) - 1.0).abs() < 1e-10);
        assert!(fast_expn(30.0) < 1e-10);
        let at_max = fast_expn(EXPN_MAX);
        assert!(at_max == 0.0, "boundary x==25 returned {}", at_max);
        for i in 0..100 {
            let x = (i as f64) * 0.1;
            let approx = fast_expn(x);
            let exact = (-x).exp();
            assert!((approx - exact).abs() < 5e-4, "x={} approx={} exact={}", x, approx, exact);
        }
    }
}
