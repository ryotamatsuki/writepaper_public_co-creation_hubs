# Stage 7 v2.1 — Saturated Partner-Surplus Accounting Correction

Date: 2026-09-10 JST

Trigger: Astra Stage-11 hostile review.

## Primitive authority

For support route `h`, gross participation benefit is `r_h` and support cost is `c~U[0,1]`. A support-side participant joins iff `c<=r_h`, subject to unit mass. Therefore

`m_h = clip(r_h,0,1)`.

Gross benefit `r_h` is not clipped. National support surplus is

`S_h = integral_0^{m_h} (r_h-c) dc = r_h*m_h - m_h^2/2`.

Equivalently:

- `r_h<=0`: `S_h=0`;
- `0<r_h<1`: `S_h=r_h^2/2`;
- `r_h>=1`: `S_h=r_h-1/2`.

Each regional government receives one half of national support surplus under the maintained symmetric regional attribution.

## Correction

The historical Stage-7 interior formula `m_h^2/4` for regional surplus is valid only when `0<r_h<1`, because then `m_h=r_h`. It is not a valid all-domain formula after participation saturates.

The corrected general regional welfare object is

`W_i = PS_i + (1/2) sum_h [r_h*m_h - m_h^2/2] - gamma*x_i^2/2`.

This formula is continuous at both clipping boundaries. At `r_h=0`, national surplus is zero on both sides; at `r_h=1`, the interior and saturated expressions both equal `1/2` nationally and `1/4` regionally.

## Effect on claims

The reported G and B3 stationary candidates are support-interior, so the corrected general formula coincides with the historical interior shortcut at those points and does not alter the local T1/T2 algebra.

The correction is material for large off-path deviations, including `x_i=0.9` and `x_i=1`, where local-public support participation saturates. Those histories must use the general primitive integral.

W1, exact private-fee transfer cancellation, is unchanged.

W2 is a local directional numerical statement at the reported G stationary candidate, with own, rival, private-profit, and direct national derivatives generated explicitly after the correction.

W3 is a numerical ranking at the reported G/B3 state pair only. It is not qualified as a ranking of globally certified or unique equilibria after the Stage-11R evidence repair.

No primitive, support cap, gross-benefit function, regional attribution rule, or welfare objective is changed. This record corrects the implementation and the domain of validity of the prior shortcut.
