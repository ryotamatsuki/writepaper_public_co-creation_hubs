# Stage 11R P1f — frozen-model primitive conditions, pointwise robustness, and novelty gate
2026-10-28 JST. Baseline is Draft PR #61 (head `cb92ec3cc5bb6fb94b20420bc5085605064d00cc`); nested exploration PR #62. **No model/manuscript/Lean/primitives are changed.** This is not a journal-ready novelty certificate.

## Executive verdict
- **Proof added:** in the frozen model, for *weak bilateral matching* the strategic cross-effects of the two governments have an explicit **primitive-level local threshold** in the regional commercialisation incidence `M_i`. Different signs produce a no-pure-strategy, unique-mixed-entry rectangle by the standard 2×2-game theorem.
- **Pointwise stronger test:** nine distinct values of `(lambda, omega_B)` have individually certified global P/P posted-fee Nash continuation, strict A strategic complementarity, strict B substitutability, and therefore a strictly positive area of public fixed costs with no pure government-entry Nash equilibrium.
- **NOT proved:** one single certified continuous rectangle in `(lambda,omega_B)`, a broad stable common government-cost rectangle, the theorem's necessity for all parameter values, or novelty over all closest predecessors.
- **Publication recommendation:** keep this material exploratory; do not replace PR #61's main contribution yet. Scholarly novelty **HOLD**.

## 1. General analytical sufficient condition, beyond a numerically fitted threshold
Let `gamma=rho*mu` tend to zero from above, **holding rho fixed** so the basic service utility `u_i(s_j)=v_i+rho*w_i` remains unchanged for all `j=P,B`. Hold the product-market primitives, `omega_i`, `lambda`, `R0`, `p_b`, `Pi`, local premium valuation `d_i=Delta_i` fixed. Assume `0<u_i<d_i`, uniform skills aligned with private taste, and the interior Cournot/success-feasibility inequalities. Write
```
ell_i=(d_i-u_i)/(2*d_i) in (0,1/2),
M_i = consumer+incumbent welfare change per successful venture,
K_i=M_i+omega_i*Pi,  S_i0=p_b*M_i+omega_i*R0.
```
The `PP` price-continuation is assumed uniquely differentiable near `gamma=0` (ensured by the near-independent monopoly strict second-order conditions/implicit function theorem). Let `D_i(gamma)=W_i(B,B)+W_i(P,P)-W_i(B,P)-W_i(P,B)` be i's government strategic cross-effect after **reoptimising both firms' prices** in each entry configuration. With `j!=i`:
```
D_i(gamma) = gamma * Gamma_i + O(gamma^2),

Gamma_i/(N_i*N_j) = M_i*A_i + omega_i*B_i,

A_i = p_b*ell_i*ell_j
      + lambda*ell_i^2*ell_j^2/4 > 0,

B_i = R0*ell_j*(d_i-3*u_i)/(4*d_i)
      + lambda*Pi *[
          ell_i^2*ell_j^2/4
          + (u_i+d_i*ell_i)*u_i*ell_j^2/(4*d_i^2)
        ].
```
This is a **weak-interaction local theorem**. Define `M_i^crit=-omega_i*B_i/A_i`. If `M_A > M_A^crit` and `M_B < M_B^crit`, then for all sufficiently small positive `gamma`, `D_A>0>D_B`. Interchanging signs produces the opposite role assignment. Strict inequalities persist on a nonzero open primitive neighborhood. A nonempty *government-cost* no-pure-equilibrium region additionally requires `T_A(P)<F_GA<T_A(B)`, `T_B(B)<F_GB<T_B(P)`; the two cost bands are exactly the two strategic cross-effect magnitudes.

### Derivation, including P/P repricing
At `gamma=0`, independent private-only sellers optimally choose `ell_i=(d_i-u_i)/(2d_i)`; the same cutoff applies in PP and when the other region is fully covered. For a private seller of region i, holding the other region in j=P or B, its best-response cutoff to first order is
```
ell_i(j) = ell_i
 +gamma*N_j*omega_i*[
     -R0*m_j/(2*d_i)
     +lambda*Pi*u_i*q_j/(2*d_i^2)
  ] + O(gamma^2).
```
The difference between `j=B` and `j=P` has `m_j(B)-m_j(P)=ell_j` and `q_j(B)-q_j(P)=ell_j^2/2`. Endogenous *feedback* from an i-fee change to j's PP participation and back to i's participation requires two matching links and thus enters only at `O(gamma^2)`. It is **not** simply assumed away. Cross-differencing local real service benefits contributes
```
N_i*N_j*omega_i*(u_i+d_i*ell_i)*
[-R0*ell_j/(2*d_i)+lambda*Pi*u_i*ell_j^2/(4*d_i^2)].
```
Cross-differencing matching/commercialisation gives
```
N_i*N_j*[S_i0*ell_i*ell_j
         +lambda*K_i*ell_i^2*ell_j^2/4].
```
Combining and using `ell_i=(d_i-u_i)/(2d_i)` yields the stated `A_i,B_i`. This proof needs a justified smooth local PP price branch; the existing global price-contraction certificate at the original numerical primitive verifies a separate, stronger pointwise case. **Do not extend the asymptotic formula without qualification to gamma=0.1**.

## 2. Parameter signs and exact baseline check
For the frozen original model, the exact weak-`gamma` coefficient `Gamma_i` (coefficient of `gamma`, not the cross-effect at `gamma=0.1`) is:

| scenario | Gamma_A | Gamma_B |
|---|---:|---:|
| heterogeneous skills, lambda=0 | -0.00000452493 | -0.000666655 |
| heterogeneous skills, lambda=0.5 | +0.001142433 | -0.000511733 |

The primitive-level market-incidence critical values `M_i^crit` at `lambda=.5` are about `-0.0532853` for A and `-0.0139830` for B; actual `M_A=-0.0157974`, `M_B=-0.0307750`. Hence A lies on the complement side, B on the substitute side **in the weak-matching limit**.

Here `Gamma_i(lambda)` is affine because `lambda` enters the first-order interaction only linearly. The exact rational zeroes are: A `lambda=0.00197258000224...`, B `lambda=2.1515873354...`. The latter is outside the original examined `[0,.5]` interval. These are **weak-gamma approximations**, not the exact baseline `gamma=.1` transition (`lambda≈.00236644` numerically). Reproduced by `primitive_robustness.py`.

## 3. More robust than one calibration, but NOT broad robustness
At **nine separately fixed** values, all price-game global BR sufficient inequalities, P/P rational root-face enclosures, and A/B strict cross-effect signs pass. No fee continuation is guessed or held fixed. Values use the original technology, distributions, downstream Cournot primitives, `N_i`, `F_P`, with only existing `lambda` and venture-profit attribution `omega_B` changed.

| lambda | omega_B | A cross-effect (approx) | B cross-effect (approx) | minimum certified cost-rectangle area |
|---:|---:|---:|---:|---:|
| .4 | .18 | +.000095055 | -.000057603 | 5.4753e-9 |
| .4 | .20 | +.000091648 | -.000054077 | 4.9558e-9 |
| .4 | .22 | +.000088240 | -.000050549 | 4.4603e-9 |
| .5 | .18 | +.000118772 | -.000055230 | 6.5597e-9 |
| .5 | .20 | +.000114746 | -.000050955 | 5.8468e-9 |
| .5 | .22 | +.000110718 | -.000046682 | 5.1684e-9 |
| .6 | .18 | +.000142506 | -.000052858 | 7.5325e-9 |
| .6 | .20 | +.000137863 | -.000047837 | 6.5947e-9 |
| .6 | .22 | +.000133217 | -.000042815 | 5.7036e-9 |

For baseline `lambda=.5,omega_B=.2`, the *full theoretical* no-pure cost rectangle widths are `D_A≈0.000114746`, `-D_B≈0.000050956`, approximately **0.099% and 0.126% of the respective government-cost midpoints**. Area `≈5.85e-9` in squared cost units. This remains narrow despite certified sign persistence under ±20% variation in lambda and ±0.02 absolute changes in omega_B. Each grid point permits its **own** cost rectangle. A shared rectangle across all nine has not been shown (and should not be claimed).

The originally published baseline public costs `F_GA=.09,F_GB=.0405` **do not** lie in the mixed-equilibrium cost rectangle. The original pure-equilibrium result remains unchanged.

## 3b. **Uniform continuous parameter-rectangle certificate** (stronger than nine separate points)

An additional standard-library rational interval proof, [`uniform_box.py`](uniform_box.py), encloses **all** original-model primitives simultaneously in

```
0.4999 <= lambda <= 0.5001
0.1999 <= omega_B <= 0.2001,    omega_A=1-omega_B.
```

The program checks, for this entire 2-dimensional closed rectangle: (i) exact sufficient inequalities for global PP fee best-reply containment and Banach contraction (<1), including all nonnegative posted-fee deviations; (ii) signed all-four Poincare–Miranda price-root faces on a *common* rational cutoff box; (iii) global optimality of the remaining private price branches; (iv) strict government strategic cross-signs after joint repricing. The actual outward rational enclosures have approximate endpoints

```
D_A in [0.00009106788, 0.00013842398] > 0,
D_B in [-0.00006020123, -0.00004171075] < 0.
```

It follows that **for every parameter point in that continuous rectangle there exists a nonempty region of government fixed costs with no pure first-stage entry equilibrium**, with a unique fully mixed government entry equilibrium (both private firms enter by strict dominance). This is a **uniform-parameter theorem**, not a grid inference, although the cost rectangle need not be identical for every point in the primitive box.

The algorithm also tested larger radii `0.001`, `0.002`, `0.01`: crude interval subtraction yielded government cross-effect enclosures that straddled zero, so those larger boxes were **not certified**. These failures do not imply a real sign reversal; they are an interval dependency/precision limitation. No wider box is claimed. The original government costs still fall outside this mixed-entry window.

## 4. Scholarly novelty gate
See [closest-predecessor audit](LITERATURE_AND_NONABSORPTION.md). The critical point is that the induced asymmetric strategic interactions are endogenous to *reoptimised private prices and locally allocated downstream new-venture surplus*, but the **game-theoretic structure** is covered by established GSH and local-public-good theory. The weak-coupling `M_i^crit` formula is a model-specific sufficient condition, not evidence that broader strategic-incidence laws were unknown. No reference-full-text theorem-level non-nesting is certified; Hua et al. (2016) remains inaccessible in full.

Do not amend PR #61 to promote this to the main theorem on current evidence. Keep Stage11 novelty **HOLD**. All checks are reproducible on the exploration branch, with CI reusing the original global P/P certificate. No journal submission, PR merge or external human review.

## 5. Reproduction
Run:
```
python hotelling_entry_experimentation/production/reviews/stage11r_pp_global_certificate.py
python hotelling_entry_experimentation/production/government_game_exploration/government_game_certificate.py
python hotelling_entry_experimentation/production/government_game_exploration/primitive_robustness.py
```
The last script saves `primitive_robustness_results.json` when run; the checked source and CI logs provide the independent exact Fraction evidence. Existing original manuscript and formal proof files are untouched.
