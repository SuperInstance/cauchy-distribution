# cauchy-distribution

A Rust library implementing the **Cauchy (Lorentzian) distribution** — its PDF, CDF, quantile function, sampling, entropy, and the Hodges–Lehmann location estimator. The Cauchy distribution is the canonical example of a "pathological" distribution: it has no mean, no variance, and violates the Central Limit Theorem.

## Why It Matters

The Cauchy distribution appears throughout physics, statistics, and signal processing:

- **Spectroscopy** — natural line shapes are Lorentzian (lifetime broadening)
- **Resonance phenomena** — the amplitude response of a driven harmonic oscillator follows a Cauchy profile
- **Robust statistics** — the Cauchy distribution is the heaviest-tailed stable distribution, used to stress-test estimators
- **Cauchy–Lorentz transform** — fundamental to the theory of stable distributions (Lévy α-stable with α=1)

The distribution is also the textbook counterexample that forces students to understand *why* the CLT requires finite variance.

## How It Works

### Probability Density Function (PDF)

$$f(x; x_0, \gamma) = \frac{1}{\pi \gamma \left(1 + \left(\frac{x - x_0}{\gamma}\right)^2\right)}$$

where $x_0$ is the location (median/mode) and $\gamma > 0$ is the scale (half-width at half-maximum, HWHM).

At the peak ($x = x_0$): $f(x_0) = \frac{1}{\pi \gamma}$

### Cumulative Distribution Function (CDF)

$$F(x) = \frac{1}{2} + \frac{1}{\pi} \arctan\left(\frac{x - x_0}{\gamma}\right)$$

The CDF has an elegant closed form — no special functions required. This makes inverse-transform sampling trivial.

### Quantile Function (Inverse CDF)

$$F^{-1}(p) = x_0 + \gamma \tan\left(\pi(p - \tfrac{1}{2})\right)$$

Sampling via $U \sim \text{Uniform}(0,1)$, then $X = F^{-1}(U)$.

### Entropy

$$H = \ln(4\pi\gamma) + \frac{1}{2}$$

This is the differential entropy — maximized among all symmetric distributions with a given scale parameter.

### Why No Mean or Variance?

The characteristic function is $\phi(t) = e^{ix_0 t - \gamma |t|}$, which is not differentiable at $t = 0$. Consequently:

- $E[X]$ does not exist (the integral $\int |x| f(x) dx$ diverges)
- $\text{Var}(X)$ is undefined
- The sample mean $\bar{X}_n$ has the same distribution as $X_1$ for all $n$ (no CLT convergence)

### Hodges–Lehmann Estimator

For robust location estimation from Cauchy samples, the Hodges–Lehmann estimator is:

$$\hat{\theta}_{HL} = \text{median}\left\{\frac{X_i + X_j}{2} : i \leq j\right\}$$

This requires $O(n^2)$ pairwise means, then a median ($O(n^2 \log n)$ with sort). It has breakdown point 0.25 and is asymptotically efficient for the Cauchy location parameter.

### Big-O Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `pdf(x)` | O(1) | O(1) |
| `cdf(x)` | O(1) | O(1) |
| `quantile(p)` | O(1) | O(1) |
| `sample()` | O(1) | O(1) |
| `samples(n)` | O(n) | O(n) |
| `hodges_lehmann_estimator(samples)` | O(n² log n) | O(n²) |

## Quick Start

```rust
use cauchy_distribution::CauchyDistribution;

let cauchy = CauchyDistribution::new(0.0, 1.0); // standard Cauchy

// Density and distribution
assert!((cauchy.pdf(0.0) - 1.0 / std::f64::consts::PI).abs() < 1e-10);
assert!((cauchy.cdf(0.0) - 0.5).abs() < 1e-10);

// Quantile roundtrip
for &p in &[0.1, 0.25, 0.5, 0.75, 0.9] {
    assert!((cauchy.cdf(cauchy.quantile(p)) - p).abs() < 1e-10);
}

// Sampling
let samples = cauchy.samples(10_000);

// Robust location estimation
let estimate = cauchy_distribution::hodges_lehmann_estimator(&samples);
```

## API

| Method | Description |
|--------|-------------|
| `CauchyDistribution::new(x0, gamma)` | Create with location & scale |
| `CauchyDistribution::standard()` | x₀=0, γ=1 |
| `pdf(x) → f64` | Probability density |
| `cdf(x) → f64` | Cumulative probability |
| `quantile(p) → f64` | Inverse CDF |
| `median() → f64` | Returns x₀ |
| `iqr() → f64` | Interquartile range = 2γ |
| `entropy() → f64` | Differential entropy |
| `sample() → f64` | Single random draw |
| `samples(n) → Vec<f64>` | n random draws |
| `hodges_lehmann_estimator(&[f64]) → f64` | Robust location estimate |

## Architecture Notes

The **γ + η = C** link: the inverse-transform sampling (γ) maps uniform randomness through the quantile function, while the arctangent-based CDF (η) provides the exact distributional mapping. Together they conserve the distributional invariant C — the generated samples have exactly the Cauchy distribution, verifiable by a Kolmogorov–Smirnov test against the analytical CDF.

## References

- Cauchy, A.-L. (1853). *Sur les résultats moyens d'observations de même nature.* Comptes Rendus, 37, 198–206.
- Lorentz, H. A. (1906). *The Theory of Electrons.* (Spectral line shapes.)
- Hodges, J. L., & Lehmann, E. L. (1963). *Estimation of location based on rank tests.* Annals of Mathematical Statistics, 34(2), 598–611.
- Nolan, J. P. (2020). *Stable Distributions.* Chapter 1: The Cauchy Distribution.
- Casella, G., & Berger, R. L. (2002). *Statistical Inference,* 2nd ed., Section 3.2.

## License

MIT
