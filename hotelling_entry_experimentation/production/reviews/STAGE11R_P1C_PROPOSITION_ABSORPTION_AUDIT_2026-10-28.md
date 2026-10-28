# Stage 11R P1c — proposition-level predecessor absorption audit (adversarial)

**Date:** 2026-10-28 JST  
**Base:** Draft PR #59 at `2f7088745f9cd3f784154bb6d433c4ef6aab187d` (frozen Stage33 model).  
**Scope:** P1–P3, C1, the three-scenario ownership comparison, and empirical implications; no new economic primitives, parameters, or simulation outcomes.  
**Verdict:** **BOUNDED PRIORITY AUDIT COMPLETED; FULL THEOREM-LEVEL NON-ABSORPTION NOT CERTIFIED; OVERALL STAGE 11 HOLD.**

## 1. A serious reducibility finding (mandatory manuscript correction)

For every admissible frozen service state, founders are either (i) absent, (ii) all participating, or (iii) uniformly distributed on the upper tail `[ell_i,1]` under private-only access. It follows EXACTLY that

```text
m_i = 1 - ell_i
q_i = (1 - ell_i^2)/2
phi(m) = m - m^2/2
q_i = phi(m_i)       for all 0, P, G, B states
q_A q_B = phi(m_A) phi(m_B)
```

The `P` service also has `Z_i = N_i [u_i m_i + Delta_i phi(m_i)]`. Hence, conditional on the discrete service state, the matching-dependent part of local welfare can be written **without a second independent aggregate state**:

```text
W_i - W_i^0 - Z_i + F_Gi 1{G active} + F_Pi 1{P active}
    = C S_i^0 m_A m_B + C lambda K_i phi(m_A) phi(m_B).
```

For the pooled full-population-mean comparator, the corresponding complementarity term is `C lambda K_i m_A m_B/4`. At FIXED participation masses, the heterogeneous-minus-pooled exposure is

```text
C lambda K_i [phi(m_A) phi(m_B) - m_A m_B/4].
```

The bracket is nonnegative for `0<=m_A,m_B<=1`; its effect on a jurisdiction's welfare follows the sign of `K_i`. At full participation the bracket is zero. Reoptimised fees and provision can change the masses, so this fixed-mass sign does **not** order equilibria.

**Adversarial verdict:** the draft's assertion that the two moments cannot be represented by one participation index is false *under its own uniform truncation assumption*. There are two payoff aggregates with different functional dependence, but not two free endogenous state dimensions. A more general type distribution, arbitrary skill–taste correlation, or selection mechanism would change this reduction; those are OUTSIDE the frozen paper. The novelty must be framed as **specific nonlinear selected-matching incidence plus downstream regional and four-entry decisions**, not as introducing a fundamentally irreducible two-index platform model.

This is a priority/interpretation correction; it neither falsifies the existing P1–P3 equations nor reruns the 16-profile equilibrium certificate.

## 2. Original-paper evidence ledger

| Predecessor | Text inspected on 2026-10-28 | Specific existing result | Overlap and residual uncertainty |
|---|---|---|---|
| Hagiu (2004b), *Two-Sided Platforms: Pricing and Social Efficiency—Extensions*, RIETI 04-E-036, https://www.rieti.go.jp/jp/publications/dp/04e036.pdf | Original 30-page PDF, Section 2, especially cutoff-integrated quality `Q=integral[q_m,q_H] q h(q)dq` and developer entry, and Section 3 Proposition 3 | Vertically differentiated developers, both participation and quality-weighted aggregation, proprietary/open quality/variety welfare ambiguities; endogenous quality investment studied in a separate extension | **Strong absorption of type-selection / aggregate-quality novelty.** Their Section 2 heterogeneity and Section 3 endogenous quality investment are separate specifications; do not pretend they already derive our exact bilateral `phi(m_A)phi(m_B)` Cournot-government theorem. |
| Hua, Chen & Zhang (2016), *Competition and Coordination in Two-Tier Public Service Systems under Government Fiscal Policy*, POM 25:1430–1448, https://doi.org/10.1111/poms.12555 | Publisher abstract and metadata only; full theorem text not open-access verified | Zero-price public provider alongside fee-charging private provider; Nash equilibrium and existence of a two-tier system; fiscal coordination | **Absorbs coexistence as a concept, and defeats claims of new free-public/paid-private competition.** Cannot determine exact necessary/sufficient conditions or exclusions from abstract; full-text theorem priority is **OPEN**. |
| Long & Wong (2009), *A Tale of Two Ports*, Review of International Economics 17:261–279, https://doi.org/10.1111/j.1467-9396.2009.00822.x | Open complete article, Sections 2–3 and 5–7 | Two geographically separated governments choose infrastructure support for local private distributors competing in service prices; private price equilibrium and strategic regional rivalry | **Absorbs the mere two-government/two-firm sequential architecture.** Their infrastructure lowers distributors' costs instead of supplying a competing zero-price basic option; they do not build founder-pair R&D success and region-specific Cournot product-entry payoffs. This difference is descriptive, not a non-embedding proof. |
| Peitz & Sato (2026), *Asymmetric Platform Oligopoly*, RAND J Econ 57:60–77, https://doi.org/10.1111/1756-2171.70008 | Open complete article, Sections 2–6; main results, including online-appendix summaries on partial coverage | Asymmetric two-group platform fees, within/cross-group network effects, exogenous and endogenous platform entry, one-sided pricing, free outside option and partial participation | **Absorbs price/network/entry/partial-coverage architecture.** Specific maintained logit/log-network specification does not literally include type-selected bilateral venture profits or four independent public service investment bits; nevertheless the `phi(m)` reduction raises a serious generalised-network interpretation. No formal non-nesting theorem established. |
| Cremer, Marchand & Pestieau (1997), *The optimal level of local public goods: A two-region model*, JPubE | Published reference already audited in Stage 11 | Local vs joint-surplus provision wedge arises from neighbour's welfare increment | **Absorbs generic logic of P3.** Endogenous matched-venture generation of the foreign term is the only model-specific part. |

**Additional caution:** Laine & Ma (2017), Kim (2024), Berliant & Fujita (2008/2011), and Lee (2026) remain material comparators, as described in the earlier Stage 11 audit. This four-paper focused audit does not make them disappear.

## 3. Proposition-by-proposition collision test

**P1: local versus aggregate marginal complementarity.** The derivative `dW_i/dlambda=C K_i/4` under full participation is a direct linear welfare-accounting identity. A local sign differing from the sum is a textbook distribution/incidence possibility, familiar from interregional spillovers and product-market rivalry. **Status: generic mechanism ABSORBED; exact application model-specific, no independent novelty theorem.**

**P2: fee, mass and mean-type selection.** With fixed partner moments and `0<U<T`, monopoly linear-demand pricing gives `f*=(U+T)/2`, `ell*=1/2-U/(2T)`. Increasing `T` increases the cutoff, decreases `m,q`, and increases the conditional mean. This is a uniform upper-tail and monopoly-pricing calculation, against a much earlier platform/quality-selection literature (especially Hagiu). **Status: generic mechanism ABSORBED; sensitivity to this model's `lambda` is a specialisation, not a priority claim.**

**C1: selected relative to random admission.** `q-m/2 = m(1-m)/2 >= 0` is basic truncation algebra. The fixed-mass comparator is NOT an equilibrium alternative selection rule implemented by the frozen game. **Status: algebraic identity, not novel.**

**P3: public addition and cross-border wedge.** `F_national-F_local=Delta W_j` is a welfare-accounting externality, not a new local-public-finance theorem. Hua, Long–Wong and Cremer et al. substantially precede its general components. **Status: generic threshold logic ABSORBED; the particular `lambda K_j ell_i^2/4` matched-project term is model-specific.**

**Finite three-scenario provider-configuration contrast.** The certified three outcomes `(B,P)->(B,B)` with skills as lambda goes 0 to .5, versus `(B,P)` under the pooled comparator at .5, hold for the documented narrow `F_GB` open interval (about 0.2632% of frozen cost) and small 3D parameter rectangle. They are a **valid frozen-model numerical existence witness**, not an all-parameter theorem, not robust to Beta(2,2), and not a verified previously unproved economic phenomenon in all nearest models. **Status: model-specific composite witness NOT MATCHED to an exact predecessor in the inspected texts; non-absorption as an original theorem remains UNPROVED.**

**Baseline global P/P price verification.** The novel element for this project is proof quality and reproducibility of the chosen continuation. A baseline-specific certification is not, itself, an economically new platform-equilibrium characterization against Peitz–Sato. **Status: technical credibility contribution, not priority claim.**

**Three empirical implications.** Each is a conditional mapping/measurement design rather than a separate original theoretical discovery or a causal result. **Status: editorially useful for JICT, not independent novelty.**

## 4. Formal nesting / reducibility stress test

* **Eliminate q:** `q=phi(m)` for every admissible frozen state; **PASS, exact**.
* **Eliminate matching technology entirely into a generic two-region participation externality:** the reduced welfare term can be reproduced by an appropriately chosen nonlinear `E_i(m_A,m_B)`; **PASS at the reduced payoff-function level**, but this is not proof a named predecessor's exact microfoundations or price equilibrium are equivalent.
* **Recover every current four-entry equilibrium from Hagiu/Peitz/Long/Hua by a parameter substitution in a stated published theorem:** **NOT SHOWN**. Their specific primitives/timing/instruments differ.
* **Show no precursor can ever reproduce the ownership reversal under any extension:** **NOT SHOWN AND NOT A DEFENSIBLE NEGATIVE CLAIM**.
* **Show the skill-vs-pooled contrast exists on a parameter-open region of the frozen game:** **PASS within the narrow certificate already recorded**, not re-established by this report.
* **Use the number of institutional players or the existence of a free public tier as a novelty argument:** **FAIL** (Long–Wong, Hua).
* **Use two independent aggregates m and q as a novelty argument:** **FAIL** (single-index reduction; Hagiu's quality-weighted integration).

## 5. Required manuscript positioning and Stage gate

The defensible claimed object is **a bounded equilibrium illustration in which fee-selected ability changes the curvature of cross-region matching benefits, those benefits are distributed through separate successful-venture Cournot markets, and four independent local provision decisions consequently differ between skill-based and mean-equivalent innovation technologies**.

Do not claim: new price screening; first quality-weighted participation; a general theorem that public provision improves welfare; two independent participation-state dimensions; a universally valid physical hub policy; general SPNE pricing uniqueness; or a robust ownership reversal. Explicitly acknowledge that even the full conjunction may be reproducible by a sufficiently general earlier framework, which is not settled by different labels or institutions.

**Stage 11R-P1c bounded audit: COMPLETE with material correction. Global novelty clearance: HOLD.** The next actual gate is an external-economist adversarial assessment of whether the narrow calibrated illustration merits a full applied-theory paper, with Hua full text / appendix-level comparison where lawfully accessible. If an earlier proposition is established to encompass the exact combined equilibrium claim, return to Stage 6 rather than adding an ad hoc primitive.

## 6. Audit provenance and independence limits

Sources are linked above; Hagiu original PDF and Wiley complete open texts (Long–Wong, Peitz–Sato) inspected directly. Hua exact proposition comparison is restricted to the publisher abstract, therefore NOT certified. No external human referee was consulted. The symbolic reduction is elementary algebra independently checked against the frozen forms and by standard-library exact-rational test cases. This report is not a comprehensive new systematic literature review, publication-priority guarantee, or a claim that all original numerical cases were rerun again.
