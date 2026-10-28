# Stage 8 — Canonical Theory Freeze (v1)

**Theory object ID:** `STAGE33-2026-10-28-THEORY-FREEZE-v1`  
**Status:** `THEORY FROZEN — GO TO REPRODUCIBILITY SETUP`  
**Gate closure:** Stage 4A has a scope-limited independent mathematical GO; Stage 6 and Stage 7.5 bounded integrated contribution GO; Stage 7.5A has `GO — GENERALITY / QUANTIFIER / PORTABILITY CERTIFICATION PASS` for the explicitly restricted claim register; Lean `FORMAL VERIFICATION PASS` covers the proof-critical core. The author explicitly self-attested completion of the three named proof/equilibrium review issues on 2026-10-28 JST. This records an author-reported scope review, not independent author-run Lean compilation or a complete mathematical reproof; see the immutable author accountability record and CI evidence below.  
Date: 2026-10-28 (JST). Workflow: `ryotamatsuki/research-paper-workflow` v2.9.  
Scope: Stage 33; *not* the separate `public_two_sided_platform_stage8_theory_freeze` model already on this repository's main branch.

## 1. Research question and positioning (locked)

**Working title:** *Public and Private Innovation Intermediation: Firm Sorting, Regional Market Effects, and Endogenous Provision*.

When firms with complementary R&D and commercialization skills in two regions can access public free basic innovation intermediation and private fee-based specialist services, under what conditions do local governments and profit-maximizing intermediaries enter or coexist, and when do decentralized provision choices diverge from nationally preferred provision?

Contribution class: **bounded integrated/model-specific economic characterization**, coupling endogenous skill composition (`m_i,q_i`), matched innovation, downstream differentiated-product competition and local-vs-national provision incentives. Not a priority claim for discovering price sorting, welfare externality, or ownership reversal for the first time. Do not assert that any complete class of earlier theorems has been ruled out.

## 2. Players, information and chronology (locked)

Two regions `i in {A,B}`, each with a local welfare-maximizing government `G_i` and a local profit-maximizing privately financed service intermediary `P_i`; total four first-stage strategic agents. Each simultaneously and independently chooses a binary provision/entry action, giving sixteen profiles `s=(s_A,s_B)`, where `s_i in {0,G,P,B}` and `B=Both (G+P)`. A government's payoff is its **own region's welfare**, not consolidated national welfare; a private operator's payoff is its fee revenue minus its own real fixed cost.

After entry, present private operators simultaneously choose posted nonnegative fees `f_i in [0,+infinity)`; founders privately know their types and select a service or outside option; independently formed A–B matched founder pairs invest in joint-venture R&D; successful ventures compete against pre-existing producers in differentiated Cournot product markets in both regions. Investment and Cournot choices are interior only on the explicitly admissible baseline domain; other branches require separate continuation analysis. Matched founder pairs have their own separate product niches, so aggregate projects do not double count a single venture's product profit.

Firm types `theta_i in [0,1]` are independent across regions and uniform in the baseline; technical ability and willingness to pay for private premium service coincide `z_i=theta_i`. This is a **substantive positive-correlation assumption**. Number of potential founders per region: `N_i>0`. Private fees are nonnegative; firm outside option is zero.

## 3. Product market, R&D and regional social accounts (locked)

With `n_i` incumbent firms and a successful entrant of quantity `y_i`, region-i representative consumer utility is `a_i Q_i + a_i^v y_i -(b_i Q_i^2+2 chi_i Q_i y_i+d_i y_i^2)/2`, using parameters ensuring a positive-definite demand Hessian and positive interior equilibrium quantities. Incumbent per-firm quantity `x_i` and entrant `y_i` are the Stage30 differentiated Cournot solutions; zero-hub incumbents operate normally. Let `Pi=Σ_i d_i y_i^2 > 0` be the JV's successful-product operating profits across both markets. Let `M_i` be the change in consumer surplus **plus incumbent producer surplus** from a successful new good in region i, excluding JV profit. Let `omega_A+omega_B=1`, with each share nonnegative.

Per-match technological success `p(e,theta_A,theta_B)=p0+zeta e+lambda theta_A theta_B`, with all relevant success probabilities in `[0,1]`; R&D cost `k e^2/2` with `k>0`. Optimal project investment `e*=zeta Pi/k` is type- and lambda-independent in this separable specification. Define `p_b=p0+zeta^2 Pi/k`, `R0=p0 Pi+zeta^2 Pi^2/(2k)`, `K_i=M_i+omega_i Pi`, `S_i^0=p_b M_i+omega_i R0`, `S_i^1=lambda K_i`.

`C=rho mu N_A N_B>0`; participation mass `m_i=Pr(join)`, selected first skill moment `q_i=E[theta_i 1{join}]`. Match-count channel: `L0=C m_A m_B`; skill-weighted channel: `L1=C q_A q_B`.

## 4. User choices, prices, active sets and transfers (locked)

Conditional basic participation benefit `u_i=v_i+rho w_i 1{region j has any hub}>0`. Define `a_i^m=rho mu N_j omega_i R0`, `b_i^m=rho mu N_j omega_i lambda Pi`; these are **matching coefficients**, not Cournot demand slopes.

No hub: `m_i=q_i=0`. Free government `G` and joint `Both`: `m_i=1,q_i=1/2`. P-only type `theta` compares `U_P=u_i+a_i^m m_j+(Delta_i+b_i^m q_j)theta-f_i` with zero. Its cutoff:
`ell_i=clip[(f_i-u_i-a_i^m m_j)/(Delta_i+b_i^m q_j),0,1]`,
`m_i=1-ell_i`, `q_i=(1-ell_i^2)/2`.

If other-region participation is fixed, put `U_i=u_i+a_i^m m_j`, `T_i=Delta_i+b_i^m q_j`. For `U_i>0,T_i>0`, globally optimal monopoly price on `f_i>=0` is `(U_i+T_i)/2` if `U_i<T_i` (interior sorting), and `U_i` if `U_i>=T_i` (full-coverage price kink). Zero demand is dominated by feasible positive-revenue prices.

If Both: all founders use free basic intermediation and premium P is bought iff `Delta_i theta_i>=f_i`. The private premium monopoly price is `Delta_i/2`, premium participant share `1/2`; direct regional welfare service benefit `N_i(u_i+3 Delta_i/8)`, private gross revenue `N_i Delta_i/4`. If P-only the direct benefit is `N_i(u_i m_i+Delta_i q_i)`.

If (P,P), prices are **simultaneously strategic**, not two isolated monopolies. Participant cutoffs solve the coupled fixed point. A global pure-price equilibrium is unique only when the **documented Stage4/4A sufficient conditions** hold; do not assert this for the entire primitive domain. For all first-stage unilateral deviations, downstream continuation must be solved or explicitly classified `UNRESOLVED`.

Private fee payments are transfers between founders and local private operators; exclude fees from aggregate real service benefit and do not subtract them twice. Real fixed costs `FG_i` and `FP_i` enter local welfare and private operators' profits as appropriate.

## 5. Welfare, planner benchmark, and first-stage equilibrium (locked)

`W_i(s)=W_i^0+Z_i(s)-FG_i 1{G_i enters}-FP_i 1{P_i enters}+C m_A(s)m_B(s) S_i^0+C q_A(s)q_B(s) S_i^1`, where `W_i^0` is no-venture incumbent consumer+producer surplus and `Z_i` is gross direct service benefit (0/G/P/Both as above).

Private payoff `pi_i^P(s)=N_i f_i(s) D_i(s)-FP_i 1{P_i enters}`, where premium demand `D_i=m_i` if P-only and `D_i=1/2` if Both; absent private provider payoff zero. A first-stage pure Nash equilibrium is a profile with all four binary unilateral action flips weakly unprofitable, evaluated using valid contingent downstream equilibria. It is not declared globally unique at indifference points. Off-path continuation classes include no hub, only free, P-only interior, P-only full coverage, Both premium, P/P coupled pricing, and any uncovered boundary or multiplicity.

The planner chooses only `s in {0,G,P,Both}^2` to maximize `W_A(s)+W_B(s)`, **subject to the downstream decentralized equilibrium and its selection**. This is `restricted-instrument provision optimum` (not `first best`). A comparison between two fixed configurations is **not** by itself the full sixteen-profile planner optimum.

## 6. Claim register and maximum defensible quantifiers (locked)

| ID | Admissible statement | Proof maturity | Formal coverage |
|---|---|---|---|
| P1 | With full participation, fixed interior differentiated Cournot and independent uniform types, `dW_i/dlambda=C K_i/4` and `d(W_A+W_B)/dlambda=C(M_A+M_B+Pi)/4`. Sign may differ across regions/nation. | `PROVED CONDITIONAL` | only auxiliary order/algebra; full market model not Lean-encoded |
| P2 | With other-region participation fixed, `0<U<T`, `theta=z`, positive lambda sensitivity of T: optimum fee and cutoff rise, participant mass and total skill moment fall and selected average ability rises. The inverse-hazard generalization is analytic with stated conditions, not a full-game universal claim. | `PROVED CONDITIONAL` | `private_fee_increases`, `cutoff_increases`, `participation_decreases`, `aggregate_skill_decreases`, skill moments: uniform special case only |
| P3 | Conditional P-only to Both addition, other region fully covered: `F_i^local=H_i+C(S_i^0 ell_i+lambda K_i ell_i^2/4)`, `F_i^national-F_i^local=C ell_i (S_j^0+lambda K_j ell_i/4)`; either sign possible. | `PROVED CONDITIONAL` | `planner_local_wedge`, `positive_foreign_externality`, `negative_foreign_externality` |
| P4 | Baseline uniform/Cournot parameter example: the conditional local public-entry surplus crosses zero once on `0<=lambda<=0.5` around `0.3511840318401`. Full 4-player profile switch (Both,P) -> (Both,Both) is asserted only for independently checked baseline continuation and strict-deviation cases; at root both may be NE. | `CONDITIONAL + NUMERICALLY SUPPORTED EXISTENCE` | `public_entry_polynomial_identity`, `public_entry_iff_polynomial`, `baseline_phi_positive`, endpoint signs, binary-flip logic; global 16-profile SPNE not Lean-proved |
| L1–L2 | General inverse-hazard sorting and free-baseline premium price, under monotone inverse hazard, regular density, interior optima, positive valuations. | `PROVED ANALYTICALLY CONDITIONAL` | not generalized in Lean |
| L3 | Continuous-type regional vs national welfare threshold difference conditional on other region full participation, documented means and moments. | `PROVED ANALYTICALLY CONDITIONAL` | uniform specialization of welfare algebra in Lean |
| C1 | Fixed participation cutoff: selection raises skill-weighted matches relative to random selection of equal mass; not a statement comparing **reoptimized equilibria**. | `PROVED CONDITIONAL` | per-region moments in Lean, two-region product identity not itself encoded |

**Numerical witness (not generic theorem):** `N_A=3,N_B=1,n_A=5,n_B=2,Delta_A=Delta_B=1,rho=.4,mu=.25,omega_A=.8,FG_A=.09,FG_B=.0405,FP_A=.65,FP_B=.24,p0=.05,zeta=.005,k=.0004,chi=.9`. At `lambda=0`: NE (Both,P), planner (Both,Both). At `lambda=.5`: NE and planner (Both,Both). At `lambda*=~.35118403184`, the two public-entry responses tie. The heterogeneous-vs-pooled success comparison at lambda=.5 also changes the decentralized first-stage ownership while preserving equal full-coverage innovation mean, but is a **numerical illustration**, not a general ownership theorem.

## 7. Portability, negative tests, equilibrium multiplicity (locked)

- Two pre-specified adverse environments: `Beta(2,2)` ability distribution, and differentiated Bertrand (with changed demand substitution). Neither preserves the same full 4-actor provision switch at the baseline fixed costs. These are `NUMERICAL ROBUSTNESS / NEGATIVE PORTABILITY`, not formal refutations of all variants.
- At lambda star (baseline) first-stage P and Both for region B can both be best responses. At `FP_B=.25` the P entry indifference produces another two-NE case; never write universal equilibrium uniqueness.
- `v_B=.9` P/P pricing lies outside the sufficient global pure-fee uniqueness region; selected numeric continuation `UNRESOLVED` and excluded from certified baseline.
- Private fees, transfers, market competition losses and constrained planner choices must be evaluated consistently across compared states.
- Stop rule: no untested Stage34 correlation mechanism, subsidy contract, new players, or physical-distance assumption enters this freeze without Stage3/4/4A/6 and affected Stage7.5A recertification.

## 8. Verified artifacts and source-boundary register (locked)

- Stage 33 source: `hotelling_entry_experimentation/STAGE_33_SECTORAL_COMPLEMENTARITY_SORTING_PUBLIC_PRIVATE_ENTRY.md`, research base SHA `2aae885a8781113683d061c918e7071d7bf62ba2`.
- Stage 4 tariff repair and independent Stage 4A audit: original reproducibility files from the previous research package, `STAGE4_REPAIRED_TARIFF_AND_EQUILIBRIUM.md`, `STAGE4A_REAUDIT_REPORT.md`. These are not silently replaced by the uncorrected Stage 33 research-branch text.
- Stage 6 integrated scope re-kill GO (prior first-discovery NO-GO remains for rejected claim); Stage 7 welfare; Stage 7.5 case for bounded integrated publication.
- Stage 7.5A prior CONDITIONAL GO is now closed by a separately dated scope/Lean/author-accountability closeout `STAGE75A_LEAN_FORMAL_CLOSEOUT_AND_REMAINING_GATE.md`; its negative Beta/Bertrand alternatives and quantifier ledger remain unchanged. Lean **formal core PASS** is a restricted certificate, not a model-completeness claim.
- Lean source `formal/Stage33Formal.lean` at SHA `61d6ecf2fa2fdf70428518d099f24f531e715b81`; 18 theorems, CI https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/40609975310 (success, 8656 jobs, proof-hole audit PASS), pinned Lean v4.32.1 and Mathlib `520045ab14e26149ee970e2e617ca04b09bde5d6`.
- Formal certificate `theorem_certificates/STAGE33_LEAN_FORMAL_VERIFICATION_CERTIFICATE.md` at `e6d0964422a78c70d26faea17b7b6991d360cf61`. Only `PROOF-CRITICAL CORE`; no unsupported FULL GAME claim.
- Exact rational coefficient provenance: `stage8_formal/exact_parameter_certificate.py` plus `stage8_formal/meta/stage33_exact_rational_certificate.json` in prior research artifact. Independent 2026-10-28 recomputation passed. Mapping from Cournot inputs to formal coefficients is external to Lean and must remain cited as such.
- The author's substantive research positioning and subsequent explicit **「著者自身の確認しました」** confirmation are documented in the conversation and `STAGE8_AUTHOR_INTELLECTUAL_CONTRIBUTION_RECORD.md`. Confirmation is self-attested for the three named issues; no authored complete proof reconstruction, independent external human peer review or author-run Lean compilation is asserted.

## 9. Institutional interpretation, prior art, exclusions (locked)

Public G = taxpayer-financed **free basic service access**, possibly operated by a private contractor; it does not uniquely identify direct public ownership or physical co-location. Stage 33 does not compare subsidy/contracting or digital-only matching technologies, estimate causal program effects, endogenize government budget finance, or guarantee a particular site is needed.

Closest work includes Hagiu (2004a/2004b), Berliant–Fujita (2008/2011), Takahashi (2004), Laine–Ma (2017), Kim (2024), Hua–Chen–Zhang (2016), Peitz–Sato (2026). The priced/open platform and sorting effects are not asserted to be new discoveries; some closest prior theorem texts remain unverified, so universal non-absorption is not claimed.

## 10. Formal Stage 8 gate closure and change control

**Final workflow verdict: `THEORY FROZEN — GO TO REPRODUCIBILITY SETUP`** for the exact bounded Stage 33 theory specified above.

### Admissible evidence for closure

1. Stage 4A: independently reconstructed payoffs, pricing and 16-profile deviations; `GO — MATHEMATICAL ADVERSARIAL CERTIFICATION PASS` only for the explicitly admitted baseline and sufficient-condition tariff domains, not the full unrestricted parameter space. Retain original Stage 4A assumptions and corner/multiplicity exclusions.
2. Stage 6: `GO` for the integrated, economically differentiated but **model-specific** contribution. The disallowed claim that the generic sorting-driven public entry reversal is a first discovery remains rejected.
3. Stage 7: restricted-instrument national welfare comparator, full 16-profile optimum only when the relevant downstream continuation is valid, and separate regional spillover conditions.
4. Stage 7.5: contribution-completeness `GO` for a diagnosis-oriented theory paper, not an implementation prescription.
5. Stage 7.5A: `GO — GENERALITY / QUANTIFIER / PORTABILITY CERTIFICATION PASS` **within the exact claim register**, after retaining adverse Beta/Bertrand portability findings and explicit equilibrium-selection conditions.
6. Formal: `FORMAL VERIFICATION PASS — PROOF-CRITICAL CORE` (Lean 4 v4.32.1, pinned Mathlib SHA `520045ab14e26149ee970e2e617ca04b09bde5d6`; 18 kernel-accepted theorems; CI run `40609975310`; no `sorry`, `admit` or project-specific axiom; `#print axioms` standard axioms only).
7. Author accountability: author explicitly reported completion of three named proof/equilibrium review points in conversation on 2026-10-28 JST, recorded as **author self-attestation** in `STAGE8_AUTHOR_INTELLECTUAL_CONTRIBUTION_RECORD.md`. This is a personal review statement, not an AI-produced approval, an external peer review, or evidence the author personally ran kernel validation.

### Freeze integrity and exclusions

The previously prepared `STAGE8_CANONICAL_THEORY_SCOPE_LOCK_RC1.md` remains immutable as a **pre-freeze candidate**. This `v1` document is the frozen canonical manuscript authority. The frozen theory is tied to the formal source at commit `61d6ecf2fa2fdf70428518d099f24f531e715b81`, its formal certificate on the same working branch, and the pinned successful CI run, with the post-freeze **integration commit SHA** to be recorded in the separate freeze manifest after research-branch merge.

No claim of complete Lean certification of Stage33 global SPNE, unrestricted P/P price equilibrium, universal welfare results, universal technical/taste correlation, physical-facility indispensability, subsidy superiority, or general first-ever reverse entry mechanism is licensed.

### Theory change control

No material change to players, timing, utility/demand, information, strategy domain, participation/fixed-point boundary, game equilibrium, planner feasible set, proposition domain, theorem quantifier, formal theorem map or contribution claim without a written change request identifying the diff and rerunning the earliest affected gates: Stage 4/4A (equilibrium or strategic change), Stage 6 (contribution/novelty), Stage 7 (welfare), Stage 7.5A (scope/portability) and the Formal Verification Gate (formal statement changes). Subsequent drafting is limited to exposition and faithful restatement. This freeze does not automatically certify new Stage34 models or results.
