# Stage 13R-TP — Theory Integration Audit

Date: 2026-09-03

## Startup authority map

- startup `origin/main`: `6f1d36d22a763b2615fcad0c9b8428c84a7651ca`
- startup open PR count: `0`
- production manuscript: `paper/main.tex` plus `sections/*.tex`
- economic-model authority: Stage 8 freeze, payload `53405fa904dd31817639a734de2063158ec69321`
- manuscript/referee authority before this stage: merged Stage 11 PR #37
- journal-fit authority: merged Stage 12 PR #38
- current analytic proof-status authority: merged L3-0 PR #39 plus L3-1 PR #40
- exact canonical verifier: `analytic_level3_reopen/code/verify_l31.py`

The Stage-12 description of the headline proof ceiling as computational/local is a historical snapshot. It remains valid for the manuscript as it existed when Stage 12 was run, but it is superseded for current proof status by L3-0/L3-1. Journal selection remains a Stage-12 authority.

## Section-by-section pre-integration audit

### Abstract — MAJOR revision required

The pre-integration abstract described the headline as a numerically certified existence statement with local open-set persistence. That wording materially understated the L3-0 analytic result and made the paper look calibration-driven. It also did not distinguish ordinary numerical verification from the L3-1 exact algebraic certificate.

Action: promote the small-beta theorem; retain explicit global-theorem limitation; describe the canonical certificate as exact computer-assisted verification rather than as the source of generality.

### Introduction — MAJOR revision required

The economic mechanism was already clear, especially the fixed-price network force and the offsetting private-price response. The weakness was proof ordering: the numerical witness appeared before an analytic theorem because the manuscript predated L3-0/L3-1.

Action: reorder the contribution as mechanism -> analytic local theorem -> exact canonical certificate -> numerical magnitudes/robustness. Preserve the early concession that endogenous follower-price effects are known components.

### Model — NO substantive change

The model already cleanly defines the two public investors, shared private intermediary, timing, project single-homing, partner multi-homing, regional welfare objectives and private profit objective. No new primitive, player, timing, objective or side is needed for the proof upgrade.

Action: leave the frozen economic structure unchanged.

### Equilibrium and benchmarks — MINOR/NO change

The B3/G distinction was already disciplined: B3 fixes `bar p_T=p_T^G` as a matched-price diagnostic and is not interpreted as price regulation. The best-response slope identity already makes the cross derivative the relevant sign object under a strict SOC.

Action: use this existing notation as the bridge to the theorem; no benchmark redefinition.

### Main results — MAJOR revision required

The old Proposition was explicitly computational. The section contained the mechanism but lacked a referee-visible theorem hierarchy.

Action implemented:

1. Proposition — first-order B3 complementarity around beta zero.
2. Proposition — exact repricing decomposition/dominance identity.
3. Theorem — local strategic sign reversal for sufficiently small positive beta.
4. Corollary — endogenous commercial adaptation can change the sign of public strategic interaction.
5. Proposition — exact certification at the canonical rational parameter vector.

### Welfare — NO substantive change

The welfare section already keeps the strategic-sign reversal distinct from the welfare externality and treats the private access fee as an aggregate transfer. No new planner or policy theorem is justified by the proof upgrade.

### Robustness — MODERATE revision required

The old section made numerical perturbations and IFT continuity do too much rhetorical work. The new architecture separates:

- analytic beta-neighborhood theorem;
- exact canonical algebraic sign verification;
- numerical perturbations for magnitude/knife-edge exposition;
- local smooth distribution/network/asymmetry robustness.

No broad parameter-plane heatmap is added. The implicit epsilon and the exact beta=0.05 certificate do not by themselves prove one connected regular branch for every beta between them, so the manuscript does not make that inference.

### Institutional interpretation — NO substantive change

The analogue-class interpretation remains appropriate. No causal private repricing evidence is claimed.

### Related literature — POSITIONING RETAINED

The existing section already concedes the closest component overlaps, especially Lopez & Vives on follower responses, and defines the surviving contribution as the public-public sign reversal under the nested matched-price comparison. A targeted September 2026 refresh did not identify a closer model that directly supplies the same two-public-investor sign reversal. Recent platform pricing/investment work continues to absorb components rather than the whole public-public comparison.

### Discussion/Conclusion — MAJOR proof-status reconciliation required

Both sections still described the main contribution as numerical existence plus local persistence. They are revised to state the analytic local theorem and exact canonical certificate while retaining the absence of a global primitive classification.

### Appendix — MAJOR revision required

The old appendix documented a floating-point production certificate. It did not expose the L3-0 proof or L3-1 exact root/sign closure.

Action implemented:

- publication-quality small-beta proof;
- explicit beta/delta normalization;
- reason matched-price/equilibrium path derivatives vanish from the B3 first-order beta derivative;
- exact coupled rational Krawczyk isolating box;
- exact interval Hessian, price-SOC and price-response bounds;
- exact BR-slope intervals;
- clear separation of exact canonical certification from primitive-space characterization;
- machine-readable verifier path instead of printing high-degree polynomial coefficients.

## Mathematical consistency check

A potential alpha-power ambiguity was audited before integration. The analytic reduction defines `delta=alpha beta` and proves

`partial M_B3 / partial delta |0 = alpha^2 d^3/(Delta_i^3 Delta_j^2)`.

Therefore

`partial M_B3 / partial beta |0 = alpha^3 d^3/(Delta_i^3 Delta_j^2)`.

The apparent discrepancy is only the delta-to-beta normalization. The production theorem uses the beta derivative.

## Theory-drift assessment

No economic-model change is authorized or needed. The stage upgrades the status and presentation of claims already proved in the merged analytic tracks. Level-3 primitive classification remains frozen.

## Integration verdict before CI

Substantive theory integration: PASS.

Remaining gates are mechanical/reproducibility gates: manuscript lint, bibliography, LaTeX build, production tests, L3-0 identities and L3-1 exact verifier on the final branch HEAD.
