# Verification Results

Symbolic and numerical verification was executed with Python/SymPy/NumPy.

## Exact identities

- `mu_B-2mu_S = -p[(1-p)l+p(2-p)m] < 0`
- no-referral: `mu_B0-2mu_S=-p[(1-p)l+m] < 0`
- `dGamma/da=(4/9)(mu_B-2mu_S)<0`
- `dGamma/dsigma=mu_B-2mu_S<0`.

## Bernstein exact sign certificate

Referral:
- normalized polynomial degree `(5,2,2)`
- negative Bernstein coefficients: `0`
- zero coefficients: `3`
- minimum strictly positive coefficient: `7/5`.

No referral:
- degree `(3,2,2)`
- negative coefficients: `0`
- zero coefficients: `1`
- minimum strictly positive coefficient: `7/3`.

Therefore both `Gamma` and `Gamma_0` are strictly negative on the admissible interior.

## 500,000 seeded draws

- referral `Gamma`: 0 positive / 500,000 negative / 0 near-zero
- no-referral `Gamma`: 0 / 500,000 / 0
- partner-surplus-zero `Gamma`: 0 / 500,000 / 0
- no-downstream referral `Gamma`: 0 / 500,000 / 0
- no-downstream no-referral `Gamma`: 0 / 500,000 / 0
- second-hub cross-region externality: 222,325 positive / 277,675 negative / 0 near-zero

Maximum sampled referral `Gamma`: `-0.0005411469077696479`.

## Direct state enumeration

10,000 parameter vectors; single-hub, both-hub/referral and both-hub/no-referral distributions.

Maximum absolute discrepancy between explicit double-state enumeration and the moment closed form:

`1.4210854715202004e-14`.
