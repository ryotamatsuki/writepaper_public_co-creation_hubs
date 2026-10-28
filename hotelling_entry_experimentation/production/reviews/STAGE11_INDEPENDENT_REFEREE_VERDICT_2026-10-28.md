# Stage 11 — Independent referee / red-team report (2026-10-28 JST)

**Review target**: 14-page Stage 10 manuscript *Public and Private Innovation Intermediation: Firm Sorting, Regional Market Effects, and Endogenous Provision*; stacked research PR #52, Stage 10 revision at `e50ef9e80748b3a26be43ff730a43bc2f245063d`. Supporting numerical source: `hotelling_entry_experimentation/production/verification/stage33_reproduce.py`; exported 48-row CSV from Stage 10 CI.

**Theory authority**: `STAGE33-2026-10-28-THEORY-FREEZE-v1` (canonical Stage 8 blob `bafd537937a5abb40c57ac3af3b59df78ce5f702`). **No frozen model, theorem or Lean file is changed by this review**.

## Independent editorial verdict

**STAGE 11: HOLD / MAJOR REVISION REQUIRED BEFORE STAGE 12.**

- **Mathematical contradiction in the stated conditional P1–P3?** NOT FOUND by reconstructing the central symbolic bridges.
- **Reported three first-stage numerical equilibria and restricted planner outcomes?** REPRODUCED from the 48 rows using a separate four-player binary deviation and sixteen-profile maximization audit.
- **Robust, nontrivial, globally characterized ownership switch?** NO. The demonstration is a valid selected-baseline *existence illustration* but near its local cost indifference surface; several perturbations reverse it.
- **Independent proposition-level novel theorem relative to nearest published predecessors?** NOT CERTIFIED. Important additions and absorption threats remain, especially Hagiu (2004b), the 2016 mixed-two-tier system, 2024 mixed-platform duopoly, and Peitz–Sato (2026).
- **Submission-ready economics contribution?** NO. A referee could describe P1 as differentiation, P2 as a linear-demand monopoly cutoff, P3 as a familiar two-jurisdiction spillover identity, and P4 as a highly tuned counterexample. What remains is a carefully integrated, policy-relevant numerical mechanism; whether that supports a journal-length paper still requires an adversarial novelty/value response.

**Do not infer a fatal mathematical error or a demonstrated published-theorem absorption from these findings.** The appropriate decision is **HOLD**, not automatic abandonment of the research. No journal selection/Stage 12 GO is authorized.

## A. Proof reconstruction and known-model-in-disguise audit

**A1 — P1 full participation (PASS, conditional).** From `W_i = W_i^0+Z_i-F_i+Cm_Am_B S_i^0+Cq_Aq_B lambda K_i` and `m_i=1,q_i=1/2`, fixing the continuation gives `dW_i/dlambda=CK_i/4` and sum `C(M_A+M_B+Pi)/4`. No general first-stage comparative static follows. This is direct differentiation of the specified separable probability technology.

**A2 — P2 tariff and selection (PASS, conditional).** Holding opposite-region moments constant, let `U=u+a m_j`, `T=Delta+b q_j`, `0<U<T`. Revenue on the partial-coverage branch equals `N f(U+T-f)/T` and global optimum (including full-coverage price kink) is `f=(U+T)/2`. With `dT/dlambda>0`, `ell=1/2-U/(2T)` rises; `m=1-ell` and `q=(1-ell^2)/2` fall; conditional mean `(1+ell)/2` rises. **The mechanism is a textbook linear-demand monopoly cutoff and truncated uniform moment; standalone novelty effectively zero.** Joint endogenous P/P pricing cannot be inferred from this calculation.

**A3 — P3 welfare wedge (PASS, conditional but generic externality).** Changing local `P→B` with `m_j=1,q_j=1/2` adds `C ell_i` matches and `C ell_i²/4` skill exposure. Thus the local addition gain is `H_i-F_Gi+C(S_i^0 ell_i+lambda K_i ell_i²/4)`; adding the *other region's* welfare adds `C ell_i(S_j^0+lambda K_j ell_i/4)`. The planner minus local threshold **is just foreign external welfare**. In recentered form `W_i = W_i^0+Z_i-F_i+A_i m_Am_B+B_i q_Aq_B`; the algebra is a standard two-public-good/interjurisdictional externality identity with endogenous matching indices. This **reduces priority** relative to Cremer–Marchand–Pestieau (1997); no new welfare-theorem class is certified.

**A4 — Four-entrant game vs generic 2×2 participation game (PASS as architecture; MAJOR editorial novelty vulnerability).** For each fixed `s_j`, the four feasible `s_i` states form independent `G_i,P_i` entry bits. The four one-flip tests are necessary and sufficient for a pure Nash profile when every off-path continuation is assigned. However, binary normal-form tests themselves are standard and the submitted Lean theorem certifies only the logical equivalence, not the economic continuation. The nontrivial content must be the endogenous `m,q` map and proof that meaningful non-knife-edge ownership regions follow.

**A5 — R&D primitive (PASS, but an assumption designed to keep effort inert).** The optimum `e*=zeta Pi/k` is independent of type and lambda. Therefore skill complementarity operates through success probabilities only; any claim that more able founders *invest more* would be false. Matched ventures are assumed to have separate niches; scaling from many matches to additive Cournot surplus needs a transparent per-niche economy interpretation.

## B. Independent numerical red-team (external-to-production CSV recomputation)

The independent auditor parsed the 48-row exported panel and, without trusting its `is_NE`/`is_planner` columns, retested **all four pure-entry action flips for every one of sixteen profiles** in each of three scenarios. In all cases, the recomputed equilibria and restricted planner agree with the recorded flags:

| Technology / lambda | Verified pure NE | Restricted planner | Weakest NE entry margin | Planner gain vs NE |
|---|---|---|---:|---:|
| skill, 0 | (B,P) | (B,B) | +0.000074928993 (G_B no-entry) | 0.000338879255 |
| skill, 0.5 | (B,B) | (B,B) | +0.000031670102 (G_B stay) | 0 |
| pooled, 0.5 | (B,P) | (B,B) | +0.000456770767 (G_B no-entry) | 0.000936816106 |

The conditional region-B public co-entry results, and independent national-vs-local accounting, are:

| Scenario | Local B addition | Foreign A effect | National addition | Local FG_B critical |
|---|---:|---:|---:|---:|
| skill, lambda=0 | -0.000074928993 | +0.000413808248 | +0.000338879255 | 0.040425071007 |
| skill, lambda=.5 | +0.000031670102 | +0.000767656889 | +0.000799326991 | 0.040531670102 |
| pooled, lambda=.5 | -0.000456770767 | +0.001393586873 | +0.000936816106 | 0.040043229233 |

The baseline `F_GB=0.0405` is just `0.000031670102` (about **0.0782%** of the fixed cost) below the upper skill-complementary sign boundary. For the *skill versus pooled lambda=.5* equilibrium ownership contrast, cost must lie approximately in **(0.040043229233, 0.040531670102)**, a width of `0.000488440869` (about **1.206%** of the fixed cost). A `+0.0001` change in `F_GB` changes the skill-model `G_B` entry incentive's sign without changing any other primitive or the price continuation. This is a **MAJOR BUT FIXABLE numerical-importance/referee-value issue**, not proof the witness is invalid.

An independent full 16-profile cost-perturbation grid found the following approximate profile transitions (step size `0.00001`; grid boundaries not exact proofs):
- skill lambda=0: (B,B) to (B,P) around `0.04043`, planner transitions around `0.04084`;
- skill lambda=.5: (B,B) to (B,P) around `0.04054`, planner around `0.04130`;
- pooled lambda=.5: (B,B) to (B,P) around `0.04005`, planner around `0.04144`.

**Evidence status:** independent enumeration validates the *published CSV arithmetic* and scenario flags, not the provenance of each original fee/welfare value. Production CI proves the implementation executes and its internal assertions pass; it is not an independently derived full 16-subgame price-equilibrium theorem. Stage 11 must not pretend otherwise.

## C. Mandatory referee attacks and severity

| Attack | Result | Severity |
|---|---|---|
| P1–P3 proof bridge, fee kink, cutoff moments, foreign spillover | Conditional algebra reconstructed; no inconsistency found | PASS (bounded) |
| Entry and restricted planner enumeration | 48 rows × 4 binary deviations + all 16 planner profiles recomputed | PASS (table-level) |
| Closest paper / theorem absorption | Known mechanisms dominate; Hagiu 2004b overlooked; full proposition mapping incomplete | **MAJOR BUT FIXABLE** |
| Relabel/variable elimination | P3 recentered as classic two-index spillover; P2 simple fee-screening | **MAJOR BUT FIXABLE** |
| Economic size and knife-edge | Numerically valid, but key margin ≈0.078% of baseline government cost | **MAJOR BUT FIXABLE** |
| Full-domain private-private price game | Only baseline algorithm and external conditional Stage4A; unrestricted solution unproved | **MAJOR BUT FIXABLE** *if inflated*, currently bounded |
| Generality and portability | Frozen Beta(2,2)/Bertrand tests reportedly fail to reproduce same switch; no Stage10 reproducible adverse scripts | **MAJOR BUT FIXABLE** |
| Strong skill--willingness alignment | Identifies types `theta=z` exactly; no arbitrary correlation or empirical validation | **MAJOR BUT FIXABLE** for broader policy inference |
| Alternative subsidy/contract or digital matching | Excluded by model; interpretation must remain provision-of-service only | Scope PASS, policy completeness **MAJOR** |
| Baseline cost sensitivity | Cost perturbation counterexample to robust switch claim; current text disclaims universality | **MAJOR BUT FIXABLE** |
| Lean mapping | 18 hypotheses-conditioned theorems and exact rational identities; not full SPNE / Cournot or P/P pricing | PASS (truthfully bounded) |
| Full staged R&D new venture/additive multiple niches | Stated assumption; needs stronger economic interpretation | MINOR/MAJOR exposition |
| Bibliographic completeness | Most items correctly cited; add Hagiu (2004b) quality/vertical-differentiation extension | **MAJOR BUT FIXABLE** |
| PDF exposition | 14 pages; inspected pages including front, model, propositions, tables, appendix, references; readable and two central tables | PASS with standard copyedit pending |

**No new evidence of an error in the frozen theorem statements was found.** Numerical narrowness and novelty issues do **not** automatically invalidate Stage 8's restricted proof certificate; they challenge the *submission value* of the paper. If a genuine incompatibility between the fee-continuation model and the canonical game is later proven, route it back to Stage 4/4A and affected Stage 7.5A/Lean gates. If exact earlier propositions subsume the advertised contribution, reopen Stage 6 rather than rewriting the Intro to evade it.

## D. Source-verified closest prior-work positioning

1. **Hagiu (2004a), RIETI DP 04-E-035:** two-sided platform prices, access, product diversity and business stealing. This alone already undercuts the claim free/open intermediation is unambiguously efficient.
   https://www.rieti.go.jp/en/publications/summary/04120004.html
2. **Hagiu (2004b), RIETI DP 04-E-036: newly promoted closest predecessor.** The extension explicitly models vertically differentiated platform users/developers, quality investment by developers, and socially excessive variety under platform pricing. Thus the paper must not imply price-induced participation composition plus innovation quality is a completely new combination. What remains potentially different is the separate `q_Aq_B` bilateral founder ability product, *member-owned new-good profits*, two regions' Cournot incidence, and four independent public/private first-stage decisions.
   https://www.rieti.go.jp/en/publications/summary/04120005.html
3. **Peitz & Sato (2026), RAND 57:60–77:** open-access publisher full text states fee/participation equilibrium analysis, asymmetric platform entry, and partial-coverage extensions. Their entry is not the same two-local-government four-binary problem; no full theorem-to-theorem equivalence proven. Their formal participation analysis raises the bar on P/P price treatment.
   https://doi.org/10.1111/1756-2171.70008
4. **Long & Wong (2009), Review of International Economics 17:261–279:** freely accessible full article has two local governments and two private distributors with endogenous infrastructure investment followed by private price competition. Hence merely counting four actors is not new. Public infrastructure supports the distributors, whereas here free public intermediation is offered alongside private premium provision.
   https://doi.org/10.1111/j.1467-9396.2009.00822.x
5. **Hua, Chen & Zhang (2016), Production and Operations Management 25:1430–1448:** publisher abstract explicitly describes free government and fee-charging private service, two-tier equilibrium and possible coordination via subsidies. Close service-design predecessor; full paywalled theorem-by-theorem comparison remains OPEN.
   https://doi.org/10.1111/poms.12555
6. **Laine & Ma (2017), JEBO 140:336–353:** accessible publisher abstract and university manuscript show heterogeneous valuations, endogenous price/quality and multiplicity in public/private competition. It already has some mixed-outcome sign changes; distinguish our cross-region matched venture technology, not just mixed-provider non-monotonicity.
   https://doi.org/10.1016/j.jebo.2017.05.012 ; https://open.bu.edu/items/a044dca1-2aee-4c81-95da-28f9c468b948
7. **Kim (2024), Manchester School 92:281–295:** mixed two-sided platform duopoly, linear-city full coverage; prices may rise or fall after public entry depending on its quality. Our entry/partial coverage/cross-region venture structure differs, not the basic competition.
   https://doi.org/10.1111/manc.12467
8. **Berliant & Fujita (2008), IER 49:1251–1295:** knowledge heterogeneity and group sizes essential to cooperation already formalized; ours is not the first heterogeneous knowledge team model.
   https://doi.org/10.1111/j.1468-2354.2008.00512.x
9. **Cremer, Marchand & Pestieau (1997), JPubE 65:23–35:** local indivisible projects, Nash versus social provision and conditional grants; P3 is structurally a special-case foreign externality identity, even though the endogenous source of spillovers is model-specific.
   https://doi.org/10.1016/S0047-2727(97)00004-2
10. **Nishimura & Okamuro (2026), HIAS-E-159:** empirical complementarity of local government soft R&D support and local support organisations, but no identification of the four-actor game or policy counterfactual; useful to motivate, not certify.
    https://hias.hit-u.ac.jp/hias_dps/hias-e-159

Where publisher full text is not freely accessible, **abstract-level evidence is not counted as a complete theorem-priority clearance**. No article was falsely attributed a finding it has not been shown to establish.

## E. Required revisions before Stage 11 could close

**R1 — Contribution re-audit (priority P0, gate Stage6).** In a literal proposition-to-proposition table, test if Hagiu (2004b), Peitz–Sato (2026), Hua et al. (2016), Long–Wong (2009), and Berliant–Fujita (2008) can reproduce our key participation-ability effect, policy threshold and ownership illustration through a nested substitution. Separate (i) generic mechanisms already covered, (ii) specific structural coupling not represented, and (iii) formal novel theorem. Keep an unresolved flag for unavailable published texts.

**R2 — Economic nontriviality (P0, Stage11).** Add a `F_GB` interval/phase chart from the existing fixed model and quantify region width without moving the baseline to a more flattering cost. State up front that the pure ownership difference is a narrow but nonempty region. If insisting on a robust qualitative theorem, this requires **Stage8 change-control**, not extra rhetoric or an invented second calibration.

**R3 — Full-16 continuation audit (P0, Stage4/4A).** For every relevant off-path `(P,P)` continuation, save convergence, pure fee global-deviation gap, price bound, fixed-point contraction assumptions, and exact corner active-set classification. The current script's grid is not a parameter-global Nash proof. Flag unresolved regions `UNRESOLVED`, not approximate as valid equilibria.

**R4 — Adverse robustness replication (P1, Stage11).** Reproduce `Beta(2,2)` and differentiated-Bertrand tests from their original scripts and publish output including *failed* ownership switches; do not call them robustness success.

**R5 — Institutional and empirical scope (P1, Stage11/10).** Discuss the `z=theta` correlation more candidly; policy instruments and physical sites are not chosen. Explain why independent product niches allow additivity and whether 'real benefit' is cardinally meaningful. Do not claim construction necessity.

**R6 — Manuscript writing repairs (P1, Stage10/13).** Add Hagiu (2004b) to references and a concise negative-priority comparison in Related Literature; avoid presenting P1/P2/P3 as first results in economic theory. Preserve boundary assumptions and sign conventions in prose.

## Stage 11 exit decision and handoff

**STATUS = HOLD: unresolved MAJOR originality/editorial and robustness/reproducibility concerns; no verified FATAL theorem contradiction.** Stage 12 journal positioning is **not authorized** by this report. On completion of R1–R4, conduct a new *independent* Stage11 audit. If strong novelty cannot be shown, consider a sharply delimited methods/theory note or a theory + empirical calibration paper rather than asserting universal relevance. Submission itself remains blocked.

This is a **review-only artifact**. No frozen model primitives, Lean code, experimental outcomes, stage-10 manuscript or workflow laws were revised.
