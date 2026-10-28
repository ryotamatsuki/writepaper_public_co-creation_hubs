# Stage 11R P1g: frozen-model systematic equilibrium and novelty audit

**2026-10-28 JST.** Frozen base manuscript: Draft PR #61 head `cb92ec3cc5bb6fb94b20420bc5085605064d00cc`. Prior extensions: PR #62 (government 2×2 cycles) and PR #63 (weak-matching primitive conditions and exact primitive rectangles). **This report does not modify the economic game, parameters of the original examples, manuscript, original proof certificates or Lean.**

## 0. The exhaustive-scope protocol — what "no omissions" means

The companion [81-row machine-readable checklist](EQUILIBRIUM_COVERAGE_MATRIX.csv) enumerates **all logically distinct classes of continuation, entry, public welfare, perturbation, and close literature identified in this audit**. Each row records a unique ID, result/status, direct evidence, closest known theory, and a specific unresolved action. `coverage_check.py` mechanically rejects duplicates, empty cells and missing categories. Count by domain: **S=9** primitives/Cournot/R&D; **C=8** founders, participation, sorting; **P=11** fees and continuations; **E=14** four actors and public/private entry; **W=12** local/social welfare; **R=12** robustness and feasibility; **L=15** predecessor-literature checks. Total **81**.

This is a bounded research-inventory completeness standard, **not** a proof that the infinite set of possible theorems or literature is exhausted. Unsolved issues are expressly logged as **未検証/OPEN**, non-general results as **条件付き/基準**, mechanisms absent from the frozen model as **モデル外**. Do not count model-excluded axes as solved by assuming another technology. All 16 first-stage profiles and 64 actor-flips were checked in the original three baseline cases, not universally over all primitives.

## 1. Findings newly deduced without modifying the model

### Candidate A — Free public access can reduce real service value even at zero public cost

Fix the opposite region at full participation and keep an existing private intermediary in region i. Let its private-only cutoff be ell in [0,1/2] at an admissible private price, `d=Delta_i>0`, and `u=u_i>0`. The **gross direct-service gain** from introducing zero-price basic public access while maintaining the premium seller is EXACTLY

```
H_i/N_i = u*ell - d/8 + d*ell^2/2.
```

Its derivative with respect to ell is `u+d*ell>0`. It has the unique positive zero
```
ell_crit = sqrt((u/d)^2+1/4) - u/d.
```
For `0 <= ell < ell_crit`, free public access **reduces direct real service benefits**, despite weakly increasing participation. This is **premium-usage displacement**: under B the premium share is 1/2, whereas in P private all participants receive the specialist benefit. In particular, if P-only already covers everyone (`ell=0`), the entire matched population and skill moments are unchanged and for a fixed fully participating opposite region the venture term cancels. Then
```
W_i(B,s_j)-W_i(P,s_j) = -N_i*d/8 - F_Gi <0
```
for every nonnegative public fixed cost, **even F_Gi=0**. This is a formal model-internal result, not a universal claim about real public options or arbitrary P/P off-path fee games. No nonabsorption is claimed: Zhou, Huang, Hsu and Guo (2023, published online 2022) prove that public-private premium/basic mixed service can yield lower welfare than private-only provision under a different queueing model. [Primary abstract](https://pubsonline.informs.org/doi/10.1287/mnsc.2022.4424).

### Candidate B — Conditional crowd-out of private revenues becomes a general floor result

For any fee-game continuation and any nonnegative opposite matching/skill moments, if `G_i=0` but a private entrant exists, a unilateral posted fee `f=u_i(s_j)+Delta_i/2` induces participation of every `theta_i>=1/2` (even after the other region's users readjust). Private equilibrium **gross revenue** is therefore at least

```
R_i(P, s_j) >= N_i*[u_i(s_j)+Delta_i/2]/2 > N_i*Delta_i/4 = R_i(B,s_j).
```

The right-hand value holds exactly under B, including any opposite region pricing and participation, since premium upgrades depend only on `Delta_i theta_i-f_i`. The strict inequality does not use fixed opponent participation and therefore strengthens the conditional PR #61 revenue comparison. It does **not** imply the private entrant exits: if `F_Pi<N_i Delta_i/4`, private entry is strictly dominant.

This follows from the maintained nonnegative cross-match returns, same fee instrument and uniform taste distribution. The broad crowding-out notion is standard and is not a priority claim; Zennyo (2020), Zhou et al. (2023), and the public-option literature are close.

### Candidate C — Public participation access may **reduce the neighboring region's participation**

Suppose both regions initially have P-only provision, with interior independent-provider equilibrium at gamma=rho*mu=0, `0<u_i<d_i`, `0<u_j<d_j`, and other regularity/feasibility conditions; at fixed rho, let mu tend to zero. Region i switches **P -> B**, region j remains P-only, and prices **are reoptimised** at every stage-2 continuation (not treated as fixed). Let `ell_i^0=(d_i-u_i)/(2d_i)`.

A Taylor expansion of the unique smooth interior P/P Nash price continuation gives

```
ell_j(BP)-ell_j(PP)
 = gamma*N_i*omega_j*ell_i^0/(4*d_j^2)
    * [ lambda*Pi*u_j*ell_i^0 - 2*R0*d_j]
   + O(gamma^2).
```

Additionally, the **neighboring private seller's equilibrium posted fee** necessarily increases to first order:
```
f_j(BP)-f_j(PP)
 = gamma*N_i*omega_j/2
    * [R0*ell_i^0 + lambda*Pi*(ell_i^0)^2/2] + O(gamma^2)
 > 0, for sufficiently small positive gamma.
```
Thus *the fee rises for every positive underlying R0/omega_j, while the participation response changes sign.* In P/P the feedback induced by the seller's price change through region i's participation is O(gamma^2), not ignored at order gamma. The independent coupled-fee solver checks **both** derivatives numerically.

The **primitive sufficient condition** for other-region participation decline for all sufficiently small gamma>0 is
```
lambda*Pi*u_j*ell_i^0 > 2*R0*d_j,	 omega_j>0.
```
The other region's firm raises fees enough that the high-skill willingness-to-pay effect dominates the direct increase in match volume. The opposite strict inequality implies an increase in the other region's participation. This is a **local asymptotic theorem**, not a claim that aggregate national matches necessarily decline, nor a full-gamma sign theorem.

With original Cournot parameters and original `p0=1/20`, this coefficient is negative. Within the *same frozen family* (not modifying functional form), replace only the existing `p0` primitive by `1/1000`, keeping other primitives; the coefficient becomes positive and `p_b+lambda<1` remains valid. The exact Fraction identities and signs are checked in `additional_identities.py`. Separately, `cross_region_feedback_numeric.py` independently solves the two-provider P/P stationary conditions and the BP monopoly continuation at gamma=1e-3, 1e-4, 1e-5; both sign cases agree and the scaled numerical error converges toward zero. **This is a validated feasible existing-primitive possibility**, not a globally certified alternate four-player SPNE or theorem-level nonabsorption.

Close price/quality-selection predecessors: Anderson & Bedre Defolie (2025) on heterogeneous quality and bilateral platform prices; Teh et al. (2023) on fee-structure adjustment; Fainmesser & Galeotti (2016) on network-based price premia; Mekonnen (forthcoming) on public-option expansion raising the private fee. None of their available abstracts establishes the exact formula above, but the general pricing mechanism is not novel by declaration.

### Candidate D — Pure/mixed government regimes with opposite strategic interactions

The exact `D_A>0>D_B` and no-pure-cost rectangle, `D_A<0,D_B<0` and two asymmetric pure equilibria, and locally explicit `M_i^{crit}` are proved/discussed in PR #62/#63. These are **model-internal** additions. A generic 2×2 best-response cycle, strictly dominant private entry at low fixed costs, and strategic complementarities/substitutes are standard game theory. No exclusive new general equilibrium theorem has been certified.

## 2. Why these do not yet close the JICT novelty gap

1. A public option lowering private-quality use / welfare **already occurs** in Zhou et al. (2023); an especially close newer paper is Mekonnen's *Luck Out or Outpay? Competing with a Public Option*, accepted for AEJ: Microeconomics, on free public-option capacity raising private price and sometimes lowering welfare. These models differ materially (queueing/capacity, endogenous quantity constraints), so exact theorem nesting remains open, but the phenomenon is not a first discovery.
2. The P/P cross-border price/participation sign is more distinctive, because it links **two autonomous jurisdictions**, founder-quality matching, and **endogenous** intermediation fees. Yet it is a weak-interaction first-order calculation, and network-effect pricing with endogenous quality composition is extensively studied. No claim of theorem-level priority is warranted.
3. P1f's non-pure government entry and strategic sign flip persist on a mathematically nonempty parameter rectangle, but the economically relevant **government-cost window is narrow**, and the original public-cost pair stays in the previously certified pure-equilibrium regime.
4. The entire frozen model permits `q_i=m_i-m_i²/2`, so it is **exactly reducible** to a nonlinear single participation index. A claim that `(m_i,q_i)` constitute independent state variables is false in this model.
5. R&D effort `e*=zeta*Pi/k` is skill-independent. Physical buildings, quality investment, endogenous public fees, capacity constraints, joint ownership bargaining and subsidisation are **not primitive policy instruments**. They must not be smuggled into a frozen-model novelty proof.

**Net academic verdict:** Candidate C is the highest-value additional **model-internal** result; D has the strongest currently certified complete first-stage equilibrium structure; A and B are analytically transparent but close to known theory. None is yet certified as a standalone nonabsorbed JICT primary contribution. Preserve PR #61 and Stage11 HOLD.

## 3. Search/source protocol and outstanding blockers

The separate [literature audit](PRIORITY_LITERATURE_REGISTER.md) uses title, DOI/URL, original/abstract access level, absorbed mechanism, nonabsorbed question, and next validation task. **Hua, Chen & Zhang (2016)** full paper/proofs remain unavailable; publisher posts a 3.30 MB supplemental DOC but accessing a supplement alone does not verify every theorem in the paper. Also high priority are full-text proposition checks of **Zhou et al. (2023)**, **Mekonnen (accepted forthcoming)** and **Anderson & Bedre Defolie (2025)** before asserting priority over public/private-option welfare and cross-market price feedback.

Additional verified literature includes Peitz & Sato (2026), Monaco & Sabarwal (2016), Zennyo (2020), Kalamov & Staal (2023), and mixed-oligopoly/R&D precursor Zhang & Li (2013). Searches covered variants of public option, public/private freemium, innovation intermediary, quality-weighted network effects, local public good strategic substitutes, public firm entry/location, R&D matching and private crowdout. No negative-search result is treated as proof of absence.

## 4. Submission gate and reproducibility

- Frozen PR #61 manuscript remains unchanged; PR #62/63 remain separate. This P1g branch contains only exploratory proof checks, ledger and written audit.
- CI re-runs original 3 benchmark price-game and all-profile rational checks; P1e government-game certificate; P1f 9-point and uniform continuum rectangle certificates; P1g exact and numeric checks; and a separate audit-table integrity test.
- **Minimum pre-submission gate:** compare C's exact conditional formula with published propositions (especially 2023 Zhou and 2025 Mekonnen), resolve Hua full-text priority or document an explicit unavoidable limitation, and decide whether a result adds enough economic significance beyond standard selection/network pricing to change PR #61.
- **Stop rule:** if all candidate results collapse to known selection, public-option crowdout or standard 2×2 best-response mechanics, do not alter frozen model or assert certified novelty. Retain a narrower JICT submission or lower journal target only after transparent referee assessment.

## 5. Links and files
- [81-row coverage ledger](EQUILIBRIUM_COVERAGE_MATRIX.csv)
- [Ledger invariant check](coverage_check.py)
- [Additional exact identities](additional_identities.py)
- [Independent numerical fee-feedback test](cross_region_feedback_numeric.py)
- [Closest predecessor ledger](PRIORITY_LITERATURE_REGISTER.md)
- [P1f certified primitive/parameter boxes](PRIMITIVE_ROBUSTNESS_REPORT.md)
- [Previous P1e government equilibrium audit](REPORT.md)
