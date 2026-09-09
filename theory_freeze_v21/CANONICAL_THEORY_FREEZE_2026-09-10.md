# Stage 8 v2.1 — Canonical Theory Freeze

Date: 2026-09-10 JST
Workflow authority: `ryotamatsuki/research-paper-workflow` v2.1 @ `b27568172d4145a2b98825a5b66a071dbfe25f36`
Certified input SHA: `eeb48a3dd76ab6f43d5de175b12c3f374746db0f` (`decision/v21-stage75-full-theory-freeze`, Stage 7.5A final PASS)
Verdict: **THEORY FROZEN — GO TO REPRODUCIBILITY SETUP**

Supersession note: the earlier Stage 8 record at `7aa9de680141fb4778f4926ee05cd95a57385971` was branched from pre-final Stage 7.5A input `91849e2893a5599fcd64e8bbbb38724dcf4e5c67`. It is superseded by this freeze record solely to restore the correct authority chain. No theorem, primitive, equilibrium concept, welfare object, numerical witness, or claim envelope is changed here.

## 1. Research question and contribution

Two regional governments choose public innovation-intermediation investments while sharing a profit-maximizing private intermediary that subsequently chooses an access price. The paper asks whether activating that private repricing policy can change the sign of the public-public strategic relationship.

Approved contribution: with the private fee held fixed at the matched full-game on-path price, the local public strategic relation can be complementary; once the same private intermediary is allowed to reoptimize its price, the same public-public strategic object can become substitutable. The novelty claim is the matched-price sign comparison, not follower pricing, two-sidedness, public-private competition, or the generic idea that downstream reactions affect upstream incentives.

## 2. Players, timing, strategies, choices

1. Governments `i=1,2` simultaneously choose `x_i in [0,1]`.
2. In G, the shared private intermediary observes `(x_1,x_2)` and chooses `p_T` to maximize access-fee profit.
3. Project types choose among nonparticipation, own public hub, nonresident rival public hub, and the private route. Allocation is the upper envelope of route utilities.
4. Support-side participation is route-specific and may multihome; masses are bounded to `[0,1]`.

B3 is a matched-price fixed-price identification benchmark: the private scalar fee is fixed at the G on-path price while governments reoptimize. B3 is not a planner problem or price-regulation policy experiment.

## 3. Baseline primitives and functional form

Active repaired witness:

- `v=.1`
- `alpha=.5`
- `beta=.01`
- `rho=.15`
- `rho_T=.05`
- `kappa_L=.27`
- `kappa_T=.02`
- `tau=.35`
- `gamma=.825`

Baseline uses uniform project types, linear support-side feedback, quadratic public investment costs, and zero real operating cost for the private intermediary in aggregate-welfare accounting. The witness is constructive, not empirically calibrated.

At the repaired witness `kappa_L+tau=.62 > v+alpha=.60`, so nonresident rival-public use is dominated by nonparticipation for every project type. This is a sufficient dominance condition, not claimed necessary for sign reversal.

## 4. Equilibrium concepts and active numerical objects

G: decentralized sequential full game with global public best responses searched over `[0,1]`, private price reoptimized after deviations, and all project routes handled by the independent upper-envelope evaluator.

B3: matched-price fixed-price public Nash benchmark, with public best responses searched globally over `[0,1]` at the scalar G price.

Active repaired witness from Stage 4A and retained by the final Stage 7.5A certification:

- `x_G = 0.8371022382025995`
- `x_B3 = 0.8258903860237495`
- `p_G = 0.018403679612460814`
- `BR'_G ≈ -0.019870`
- `BR'_{B3} ≈ +0.003968`

The repaired witness is an all-regime computational global-equilibrium existence witness under the baseline functional form and documented tolerance. It is not a uniqueness theorem, primitive-space classification, or generic result.

## 5. Proposition and theorem register

### T1 — First-order B3 complementarity

Status: **PROVED, local sufficient-condition theorem**.

On a regular interior B3 stationary branch through `beta=0`, with the stated positive quality gaps, `d>0`, strict local SOC, nonsingular Jacobian, and smooth continuation,

`M_B3(0)=0` and `d M_B3/d beta |_{0} > 0`.

Therefore public local best-response slope is positive for sufficiently small positive beta on that continuing stationary branch.

No global equilibrium, uniqueness, explicit epsilon, or arbitrary-function-class claim is licensed by T1.

### T2 — Beta-zero full-game cross effect

Status: **PROVED, local closed-form sign result**.

At a symmetric regular beta-zero full-game stationary state with `T>0`, `Delta>0`, strict local SOCs, nonsingularity, and the baseline functional form,

`M_G(0) = -3 T alpha^2 kappa_L^2 / [16 Delta (Delta+T)^3] < 0`.

### T3 — Local strategic sign reversal

Status: **PROVED, local sufficient-condition theorem on regular stationary branches**.

There exists `epsilon>0` such that for every `beta in (0,epsilon)` on the assumed continuing regular B3 and G stationary branches,

`BR'_{B3} > 0 > BR'_G`.

T3 is not a global equilibrium theorem.

### T4 — Repaired global-equilibrium witness

Status: **NUMERICALLY / COMPUTATIONALLY CERTIFIED EXISTENCE WITNESS**.

At the repaired primitive vector, independent all-regime computation finds global public best responses in G and B3 within the documented certification tolerance, after reoptimizing downstream private pricing and participation for finite deviations. Opposite local slopes hold at those global best responses.

### T5 — Nonresident rival-public dominance

Status: **PROVED sufficient condition**.

If `kappa_L+tau > v+alpha`, the nonresident rival public route is dominated by nonparticipation for every project type and every history. Necessity is not claimed.

### W1 — Fee transfer cancellation

Status: **PROVED accounting identity**.

With zero real private operating cost, `-p_T n_T^F + Pi_T = 0` in aggregate welfare.

### W2 — Local coordination wedge

Status: **NUMERICAL BASELINE RESULT ONLY**.

At the repaired G witness, `dW^N/dx_i ≈ +0.49045`. This licenses only the phrase “local under-provision in the coordinated `+x_i` direction.” No first-best, global social optimum, optimal subsidy, or global underinvestment claim is permitted.

### W3 — G versus B3 welfare ranking

Status: **NUMERICAL BASELINE RESULT ONLY**.

At the repaired witness, `W_N^G ≈ .596788 > W_N^B3 ≈ .585928`. No general welfare dominance is claimed.

## 6. Archived / rejected evidence

The old vector `(beta=.05, gamma=.9, tau=.05)` is rejected as a global-equilibrium witness because a profitable finite public deviation exists once omitted route regimes are restored.

The old exact rational Krawczyk certificate remains valid only as an exact local stationary-root and derivative-sign diagnostic inside its isolating box. It is **not** equilibrium/SPNE authority and may not appear as the headline equilibrium certificate.

The old `±0.5%`, 20-draw robustness around that rejected vector is non-authoritative and excluded from the active evidence chain.

## 7. Approved robustness / generality scope

Approved:

- analytic local persistence in beta under the exact theorem assumptions;
- repaired one-vector all-regime computational witness;
- route dominance under the explicit sufficient inequality;
- exact transfer-cancellation identity.

Not claimed:

- arbitrary distribution robustness;
- arbitrary nonlinear-network robustness;
- heterogeneous-region equilibrium theorem;
- global primitive-space reversal theorem;
- uniqueness of the repaired equilibrium;
- necessary-and-sufficient parameter characterization;
- empirical calibration of `tau=.35` or the repaired vector.

## 8. Welfare / benchmark register

`W^N=W_1+W_2+Pi_T` is an aggregate accounting measure, not a solved planner problem.

Permitted labels: `decentralized full-game equilibrium`, `matched-price fixed-price identification benchmark`, `aggregate welfare`, `local coordination wedge`, `local under-provision in the coordinated +x_i direction`.

Prohibited labels absent a new certified planner problem: `first best`, `social optimum`, `optimal coordinated investment`, `optimal subsidy`, `welfare dominance theorem`.

## 9. Continuation / globality register

Stage 4A independent evaluator: `stage4a_v21_repaired/code/independent_repaired_audit.py`.

It does not import the Stage-4 repair solver. It reconstructs route allocation, participation, welfare, private pricing, public finite deviations, and local derivatives from primitives.

Audited history/deviation classes include public boundaries, the old dangerous low-investment region, detected maxima, repaired candidates, route switches, and private-price reoptimization. `UNRESOLVED` and `MULTIPLE_EQUILIBRIA` are fail-closed states rather than being coded as unprofitable deviations.

No unresolved continuation is authorized for the headline repaired witness.

## 10. Counterexample / regression register

Permanent counterexample: old low-friction vector permits a large public finite deviation/free-riding route and therefore invalidates the old global-equilibrium interpretation.

This counterexample is retained as workflow evidence and motivates the all-route fail-closed Stage 4A audit.

## 11. Literature positioning frozen

Closest component literatures include two-sided platform pricing/investment, mixed public-private competition, interregional public investment, sequential investment with endogenous follower responses, and strategic interaction altered by downstream feedback.

The paper does not claim those components are new. Surviving novelty is the matched-price public-public strategic sign reversal when the shared private repricing policy is activated.

## 12. Institutional interpretation frozen

Public innovation-intermediation and paid commercial innovation-community analogues support the plausibility of the architecture. The model does not establish that named real operators solve the exact one-period profit problem or that a causal public-investment-to-private-price response has been observed. The predicted negative private price response is an empirical prediction.

The repaired witness does not represent direct project poaching between public hubs because the nonresident public route is dominated under the active sufficient condition.

## 13. Claim-scope authority

Canonical Stage 7.5A scope ledger: `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md` @ input SHA `eeb48a3dd76ab6f43d5de175b12c3f374746db0f`.

Any later manuscript wording must respect that ledger. In particular, local stationary results may not be called global equilibria, and the repaired computational witness may not be upgraded to a general theorem.

## 14. Proof / numerics synchronization

The authority hierarchy after this freeze is:

1. analytic theorem statements and proofs certified through Stage 7.5A;
2. repaired independent all-regime global-equilibrium witness;
3. exact Krawczyk object only as a local stationary-root diagnostic;
4. repaired-witness welfare accounting and local wedge;
5. bounded computational stress tests within their stated support;
6. frozen literature-positioning ledger.

No rejected old witness, obsolete robustness draw, or local certificate may override a higher-ranked object.

## 15. Stage 8 delta audit

Relative to the final Stage 7.5A certified input, this Stage 8 changes no model object and adds no economic result. It only freezes and cross-indexes the already certified architecture.

Frozen parameterization, welfare object, theorem quantifiers, benchmark interpretation, numerical witness, robustness envelope, evidence hierarchy, and forbidden claims are inherited unchanged from Stage 7.5A.

The prior Stage 8 attempt is superseded because its parent predates the final Stage 7.5A merge. That is a provenance correction, not a theory revision.

## 16. Change control

Any post-freeze change to players, timing, utilities, route set, objectives, parameter restrictions, theorem statements, equilibrium/globality claims, or benchmark definitions invalidates this freeze.

- equilibrium/globality change → reopen Stage 4 and 4A;
- theorem quantifier/generality change → reopen Stage 7.5A and any earlier affected stage;
- welfare optimization/benchmark change → reopen Stage 7 or Stage 4 as appropriate;
- novelty expansion → reopen Stage 6.

Exposition-only changes and reproducibility engineering are permitted downstream only if they preserve this freeze exactly.

## 17. Final gate

All Stage 8 freeze dimensions are explicitly bounded. Uncertified generality blocks are excluded from the claim set rather than silently passed. Solver failure remains distinct from multiplicity and economic nonexistence. The repaired witness and the local theorem are not conflated.

**Final verdict: THEORY FROZEN — GO TO REPRODUCIBILITY SETUP.**
