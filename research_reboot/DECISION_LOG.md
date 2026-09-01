# Research Reboot — Decision Log

Project: `ryotamatsuki/writepaper_public_co-creation_hubs`

Date: 2026-09-01

Project starting main SHA: `cf2b1cc0510a15784f36f7d5e43ebb3ffbcd3668`

Workflow reference main SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`

Branch: `research-reboot/stage0-2-economic-foundations`

---

## Stage 0 — Idea / Motivation Intake

### Question

Is there an economic research question behind publicly supported co-creation hubs once the institution is separated from the legacy explanation?

### Evidence examined

- legacy JRS paper introduction/model/related literature;
- regional innovation systems and organizational-thinness literature;
- local-buzz/global-pipelines literature;
- innovation-intermediary literature;
- public intermediary/Japanese evidence;
- living-lab/orchestration literature;
- public economics/spillover literature;
- anti-confirmation evidence on substitution, policy failure, and unintended effects.

### Killed claims

- `hub + pipeline is automatically good`;
- unconditional complementarity of local and external collaboration;
- public provision superiority as a starting premise;
- physical hub as the economic primitive.

### Surviving question

> Under what conditions does decentralized regional provision of innovation intermediation fail to implement the socially efficient mix of local matching and extra-regional knowledge access, and does any resulting wedge depend on a genuinely intermediation-specific mechanism rather than a standard local-public-good spillover?

### Blockers

None sufficient to stop source audit.

### Canonical verdict

**GO**

Routing: `GO TO STAGE 1`.

### Next-stage contract

Audit the legacy source and mathematics only. Do not repair it or add a new mechanism.

### Stage-0 commits

- `5aebc1b710215a039b609a2f85552def33e966c0` — Stage 0 motivation intake
- `b31461bab1dc30107b125d80d6c8cff2df8451fc` — initial literature ledger

---

## Stage 1 — Source & Mathematical Audit

### Question

What does the legacy JRS model actually derive, and which headline results are assumption-driven rather than endogenous?

### Evidence examined

Legacy source at:

`jrs_hub_paper_revised_source_from_attached_folder/jrs_hub_paper/`

including model/results/extensions and bibliography, plus independent symbolic re-derivation from primitives.

### Mathematical result

No correction-grade algebraic error was found in the audited core formulas.

Verified from primitives:

- hub-alone welfare difference;
- fixed-cost and size thresholds;
- short-run/dynamic decomposition;
- pipeline optimum and SOC;
- maximized pipeline surplus;
- bundle threshold;
- key comparative statics.

### Killed / downgraded inherited claims

- hub–pipeline complementarity is imposed by `b*h`;
- pipeline-only is excluded by definition;
- optimal pipeline is positive under maintained primitives by construction;
- pipeline adds positive surplus rather than repairing the local redundancy state;
- redundancy is reduced form, not an endogenous repeated-network state;
- private/social wedge is an objective-function assumption, not a decentralized equilibrium result;
- Region B is passive;
- thickness and external-pool results are largely pair-count/scale effects;
- absorptive-capacity extension reweights surplus rather than forming capability endogenously.

### Surviving questions

Whether a distinct innovation-intermediation mechanism can create a decentralized welfare wedge beyond generic spillovers, standard search/matching, or generic capability formation.

### Blockers

None sufficient to prevent prior-art audit.

### Canonical verdict

**GO**

Routing: `GO TO STAGE 2`.

### Next-stage contract

Freeze the audited legacy representation. Search/compare literature; do not change the model to manufacture novelty.

### Stage-1 commits

- `cb0b881284e113bf7307b1e9e091d79a80151e4a` — symbolic verification script
- `5285acf7b894c4d736adcd48b4656f561141178a` — verification results
- `3f754c26c366429cdc0c9fae5ef50b076bffb864` — Stage 1 audit report

---

## Stage 2 — Literature Frontier / Novelty Kill Gate

### Question

Does any model/proposition-level distinction survive serious comparison with formal innovation-intermediary theory, public economics, regional innovation systems, and local/external knowledge-link literatures?

### Evidence examined

- 40 distinct ledgered papers;
- 20 serious closest/comparator papers in the closest-paper matrix;
- seminal to 2026 frontier coverage;
- working-paper/published-version reconciliation where identified;
- backward/forward/same-author neighborhood searches;
- formal intermediary theory including Hoppe & Ozdenoren, Macho-Stadler et al., Calcagnini et al., and Carayol & Sterzi;
- regional public-good theory including Cremer et al. and Wellisch;
- buzz/pipeline, RIS/thinness, absorptive-capacity, public intermediary incentives, and living-lab literature.

### Strongest prior-art finding

The old paper substantially underweighted formal economics of innovation intermediation.

The strongest single threat is **Calcagnini et al. (2019)** because it already gives endogenous TTO–firm matching/search in a noncooperative model, identifies externalities, and derives inefficient technology-transfer outcomes.

The strongest combined kill set is:

- Hoppe & Ozdenoren (2005): formal intermediation, matching, screening, participation, coordination failure, intermediary competition;
- Calcagnini et al. (2019): endogenous intermediary search/matching externality;
- Cremer et al. (1997): strategic jurisdictions, spillovers, Nash vs planner;
- Bathelt et al. (2004): local vs external knowledge architecture.

### Legacy claim decisions

1. Public hub corrects local coordination failure — **KILL AS MAIN CONTRIBUTION**.
2. Thin regions particularly need hubs — **KILL**.
3. Hub and pipeline are complementary — **KILL**.
4. Decentralized actors underinvest in cross-regional connectivity — **KILL AS GENERIC RESULT**.
5. Local interaction creates dynamic redundancy/lock-in — **KILL AS LEGACY THEORETICAL CONTRIBUTION**.
6. Publicly supported intermediaries improve welfare — **KILL**.
7. Hub-only / pipeline-only / bundle / neither partition — **KILL IN LEGACY FORM**.
8. Endogenous absorptive capacity as the rescue mechanism — **KILL AS NOVELTY BY ITSELF**.
9. Short-run KPI versus long-run welfare — **KEEP ONLY AS MOTIVATION / POSSIBLE ORGANIZATIONAL-MECHANISM INPUT**, not as a theorem.

### Initial Stage-2 blocker

Provisional status: `CONDITIONAL GO`.

Single blocker:

> Does existing formal work already combine innovation-intermediary search/matching with strategic regional/local-government provision, interjurisdictional spillovers, and planner comparison in one model?

### One authorized blocker-resolution search

A targeted intersection search was run and recorded in `search_queries.md`.

No exact joint model was identified among the strongest accessible records. This absence was not treated as novelty evidence. Instead, model-level comparison showed:

- closest intermediary theories do not contain the proposed strategic regional service-portfolio choice;
- closest fiscal-federalism models do not contain intermediary search/matching;
- closest local/external knowledge papers do not contain formal regional equilibrium/welfare analysis.

### Surviving research space

`POTENTIALLY NOVEL`, not verified novel:

> A regional innovation intermediary may supply economically distinct local-search and extra-regional-access services, and decentralized institutional/jurisdictional choices may distort the composition of those services. A theory contribution exists only if a new strategic interaction makes that composition inefficient in a way not reducible to standard matching externalities or generic regional spillovers.

### Canonical final verdict

**GO**

Routing: `GO TO STAGE 3 — CANDIDATE MECHANISM SEARCH`.

This GO does **not** approve the old model and does **not** certify theorem-level novelty. It authorizes only mechanism search.

### Next-stage contract

Stage 3 must compare genuinely distinct mechanisms. It may not resurrect:

- exogenous private myopia;
- quadratic redundancy as the headline primitive;
- imposed `b*h` complementarity;
- generic spillover underinvestment;
- generic intermediary search-cost reduction;
- generic endogenous absorptive capacity;
- fixed-cost thresholds;
- physical-hub essentiality.

Any Stage 3 candidate must generate a new feedback loop or theorem relative to Hoppe–Ozdenoren, Calcagnini et al., Cremer et al., and Bathelt et al. Otherwise return `NO-GO`.

### Stage-2 commits

- `6d646267d54eeab93a7c4d9c8e6a84ba4c9e7658` — search queries
- `07c2911f6b0e5a3207f83a6a9a6c6a480148633b` — closest-paper matrix
- `defb269cba259f796f3ca559b71194d42ba4af57` — Stage 2 report
- `8504594a85ef4e3ffdcabf525657b624b736cf94` — expanded/corrected 40-paper literature ledger

---

## Current project status after Stage 2

- Old manuscript: **not approved for revision/submission**.
- Legacy theory claims: **killed as contribution claims**.
- Economic motivation: **survives**.
- Theory-paper route: **conditionally justified for one Stage 3 mechanism-search round**.
- Empirical/mixed pivot: remains a valid fallback if Stage 3 fails.
- Stage 3 model construction: **not performed in Stage 2**.
