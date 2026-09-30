//! Cauchy Distribution
//!
//! Implementation of the Cauchy (Lorentzian) distribution PDF, CDF, quantile,
//! and sampling utilities.

use std::f64::consts::PI;

/// Cauchy distribution with location `x0` and scale `gamma`.
pub struct CauchyDistribution {
    x0: f64,
    gamma: f64,
}

impl CauchyDistribution {
    pub fn new(x0: f64, gamma: f64) -> Self {
        assert!(gamma > 0.0, "scale must be positive");
        Self { x0, gamma }
    }

    /// Standard Cauchy distribution (x0=0, gamma=1).
    pub fn standard() -> Self {
        Self::new(0.0, 1.0)
    }

    /// Probability density function.
    pub fn pdf(&self, x: f64) -> f64 {
        let z = (x - self.x0) / self.gamma;
        1.0 / (PI * self.gamma * (1.0 + z * z))
    }

    /// Cumulative distribution function.
    pub fn cdf(&self, x: f64) -> f64 {
        0.5 + (1.0 / PI) * ((x - self.x0) / self.gamma).atan()
    }

    /// Quantile (inverse CDF) for probability p ∈ (0, 1).
    pub fn quantile(&self, p: f64) -> f64 {
        assert!((0.0..=1.0).contains(&p), "p must be in [0, 1]");
        self.x0 + self.gamma * (PI * (p - 0.5)).tan()
    }

    /// Median (equals x0 for Cauchy).
    pub fn median(&self) -> f64 {
        self.x0
    }

    /// Interquartile range.
    pub fn iqr(&self) -> f64 {
        2.0 * self.gamma
    }

    /// Sample from the distribution using inverse transform.
    pub fn sample(&self) -> f64 {
        let u: f64 = rand::random();
        self.quantile(u)
    }

    /// Generate n samples.
    pub fn samples(&self, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.sample()).collect()
    }

    /// Entropy of the Cauchy distribution.
    pub fn entropy(&self) -> f64 {
        (4.0 * PI * self.gamma).ln() + 0.5
    }
}

/// Estimate location (x0) from samples using the Hodges-Lehmann estimator.
pub fn hodges_lehmann_estimator(samples: &[f64]) -> f64 {
    let n = samples.len();
    let mut pairwise: Vec<f64> = Vec::with_capacity(n * (n + 1) / 2);
    for i in 0..n {
        for j in i..n {
            pairwise.push((samples[i] + samples[j]) / 2.0);
        }
    }
    pairwise.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = pairwise.len() / 2;
    if pairwise.len() % 2 == 0 {
        (pairwise[mid - 1] + pairwise[mid]) / 2.0
    } else {
        pairwise[mid]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pdf_standard_at_zero() {
        let c = CauchyDistribution::standard();
        assert!((c.pdf(0.0) - 1.0 / PI).abs() < 1e-10);
    }

    #[test]
    fn test_cdf_standard() {
        let c = CauchyDistribution::standard();
        assert!((c.cdf(0.0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_quantile_roundtrip() {
        let c = CauchyDistribution::new(2.0, 3.0);
        for &p in &[0.1, 0.25, 0.5, 0.75, 0.9] {
            assert!((c.cdf(c.quantile(p)) - p).abs() < 1e-10);
        }
    }

    #[test]
    fn test_iqr() {
        let c = CauchyDistribution::new(0.0, 2.0);
        assert!((c.iqr() - 4.0).abs() < 1e-10);
    }
}

/// FNV-1a 64 — the digest every substrate in the SuperInstance fleet agrees on.
pub const FNV_OFFSET: u64 = 0xcbf29ce484222325;
pub const FNV_PRIME: u64 = 0x100000001b3;

#[inline]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = FNV_OFFSET;
    for &b in bytes {
        h = (h ^ b as u64).wrapping_mul(FNV_PRIME);
    }
    h
}

/// True if this crate's FNV-1a still agrees with the rest of the fleet.
pub fn canary_holds() -> bool {
    fnv1a64("café Δ 日本語".as_bytes()) == 0x024a555471370b18d
}
