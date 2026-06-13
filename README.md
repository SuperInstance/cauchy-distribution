# Cauchy Distribution

**A Rust library for the Cauchy (Lorentzian) distribution** — implements PDF, CDF, quantile function, sampling, and robust location estimation via the Hodges-Lehmann estimator.

## Why It Matters

The Cauchy distribution is the "pathological" probability distribution — it has **no mean and no variance** (both are undefined). It's named after Augustin-Louis Cauchy, and its bell shape resembles the Gaussian but with much heavier tails.

Why care about a distribution with no mean? Because the Cauchy distribution models **resonance** phenomena in physics: spectral line broadening in spectroscopy, resonance peaks in electrical circuits, and energy distributions in particle physics. The PDF is also called the Lorentzian function in these contexts.

In statistics, the Cauchy is the textbook example of why "average" doesn't always work: the sample mean of Cauchy random variables doesn't converge to anything (the Central Limit Theorem fails because the variance is infinite). Instead, you must use robust estimators like the median or the **Hodges-Lehmann estimator**.

Key properties:
- **Median** = location parameter x₀
- **IQR** = 2 × scale parameter γ
- **Entropy** = ln(4πγ) + ½ (maximum entropy among all distributions with the same IQR)
- **Stable distribution**: the sum of two Cauchy variables is Cauchy (like Gaussians)

## How It Works

**PDF** `f(x) = 1 / (πγ(1 + ((x−x₀)/γ)²))`: The characteristic bell shape with polynomial (not exponential) tail decay. At the peak (x = x₀), `f(x₀) = 1/(πγ)`.

**CDF** `F(x) = ½ + (1/π)·arctan((x−x₀)/γ)`: The arctangent function gives a smooth S-curve. This closed-form CDF enables efficient inverse-transform sampling: sample U ∈ [0,1], then compute `x = x₀ + γ·tan(π(U−½))`.

**Quantile function** `Q(p) = x₀ + γ·tan(π(p−½))`: The inverse CDF. At p=0.5, Q = x₀ (the median). At p=0.25, Q = x₀ − γ (first quartile). The IQR = Q(0.75) − Q(0.25) = 2γ.

**Hodges-Lehmann estimator**: A robust location estimator computed as the median of all pairwise averages `(x_i + x_j)/2` for i ≤ j. Unlike the sample mean, it has good efficiency for Cauchy data and 50% breakdown point (half the data can be outliers).

## Quick Start

```rust
use cauchy_distribution::{CauchyDistribution, hodges_lehmann_estimator};

let dist = CauchyDistribution::new(0.0, 1.0); // Standard Cauchy

// Distribution functions
println!("PDF at 0: {:.4}", dist.pdf(0.0));       // 1/π ≈ 0.3183
println!("CDF at 0: {:.4}", dist.cdf(0.0));       // 0.5
println!("Median: {}", dist.median());             // 0.0
println!("IQR: {}", dist.iqr());                   // 2.0

// Quantile (inverse CDF)
println!("Q(0.95) = {:.4}", dist.quantile(0.95));

// Random sampling
let samples = dist.samples(1000);

// Robust location estimation
let estimate = hodges_lehmann_estimator(&samples);
println!("Hodges-Lehmann estimate: {:.4}", estimate);
```

## API

- **`CauchyDistribution`** — x₀ (location), γ (scale)
  - `new(x0, gamma)`, `standard()` — Constructors
  - `pdf(x)`, `cdf(x)`, `quantile(p)` — Distribution functions
  - `median()`, `iqr()`, `entropy()` — Properties
  - `sample()`, `samples(n)` — Inverse-transform sampling
- **`hodges_lehmann_estimator(samples)` → `f64`** — Robust location estimate

## Architecture Notes

Provides the heavy-tailed distribution primitives for SuperInstance statistical analysis. Used in spectral analysis, robust estimation pipelines, and as a stress test for statistical algorithms (since the Cauchy breaks CLT-based methods). See the [architecture overview](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## License

MIT
