# Shared-Core Theory Route — Decision Log

Project: `ryotamatsuki/writepaper_public_co-creation_hubs`

Workflow: `ryotamatsuki/research-paper-workflow`

Workflow SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`

Stage-6 branch: `shared-core/stage6-post-hardening-novelty-rekill`

Stage-6 pivot base SHA: `bd4827136d38f2453b40030c67bb3d21a4fce6a9` (head of the prior Search-route PR #10; PR #8–#10 are preserved and not modified)

Date: 2026-09-02

---

## Pivot provenance

The prior Search/Matching route ended with `NO-GO — SEARCH ABSORBED`. The new shared-core route is a distinct theory pivot rather than a rescue of the killed search architecture. Its label-stripped object is a bilateral entry game in which two peripheral entrants compete against a shared residual core whose user value depends on the number of users remaining in that core.

The route was developed under the canonical workflow logic:

1. minimal bilateral-entry model;
2. proposition-level novelty re-kill;
3. one diagnosed hardening step, `beta -> (beta_M,beta_L)`, to isolate whether the network feedback must operate at the shared core or can instead operate only at peripheral entrants;
4. post-hardening Stage 6 novelty re-kill on the resulting theorem.

The branch intentionally does not overwrite the earlier reboot branches or PRs.

---

## Frozen post-hardening propositions entering canonical Stage 6

### Proposition 1 — Core-induced strategic reversal

Let `D0(A_M)` be peripheral establishment benefit when the other peripheral entrant is absent, `D1(A_M)` when it is present, and `Gamma(A_M)=D1-D0`. In the canonical core-network-only model,

`Gamma(A_M) = [Omega^2 - (A_L-A_M-mu)^2] / [2(1-2 beta)^2]`,

where

`mu = tau + 4 beta(1-beta)`

and

`Omega^2 = 2 beta [2 beta(1-2 beta)^2 + (1-beta) tau^2]`.

A nonempty open parameter set admits a unique interior threshold `A*` at which `Gamma` changes from negative to positive.

### Proposition 2 — Network-location diagnostic

In the two-beta diagnostic model,

`u_M = A_M + beta_M n_M`,

`u_i = A_L + beta_L n_i - c`,

network feedback confined to peripheral entrants (`beta_M=0`, `0<beta_L<1/2`) leaves the bilateral entry externality negative throughout the canonical interior domain. By contrast, core-only feedback (`beta_M>0`, `beta_L=0`) supports a negative-to-positive crossing on a nonempty open set.

### Proposition 3 — Entry-game implication

For intermediate fixed costs, `Gamma<0` yields asymmetric one-entry equilibria (anti-coordination), while `Gamma>0` yields the symmetric `(0,0)` / `(1,1)` coordination structure. Stage 6 treats this as an implication unless the literature audit establishes independent novelty.

---

## Canonical Stage 6 — Post-Hardening Novelty Re-Kill

### Question

Is the post-hardening theorem already known: specifically, does existing theory establish that bilateral peripheral-entry complementarity can arise from network feedback at a shared residual incumbent/outside option while entrant-side network feedback alone cannot generate it in the corresponding canonical class?

### Search discipline

The audit strips public, regional, innovation-hub, and metropolitan labels. It searches platform entry, installed-base competition, compatible versus proprietary networks, free entry under network effects, endogenous outside options, local public-facility entry, and common-rival/common-enemy theory.

### Strongest newly identified threat

Amir, Evstigneev and Gama (2021), *Economic Theory*, compares firm-specific proprietary networks with a single compatible industry network. Their Proposition 3 establishes opposite entry/viability effects: starting from monopoly, additional firms improve viability in the single-network model but reduce viability in the firm-specific-network model.

This kills any broad contribution claim of the form `where/how network effects operate changes the effect of entry`. That is prior art.

It does not, however, contain the post-hardening bilateral residual-core theorem. In their single-network case, additional firms jointly enlarge the same compatible network. In the candidate model, entry by a peripheral platform removes users from a residual incumbent core, lowering a shared outside-option value faced by the other peripheral entrant. Their object is industry viability as the number of firms changes, not a named entrant-to-entrant cross-entry payoff `Gamma`.

### Other decisive threats

- Peitz and Sato (2026): asymmetric platform oligopoly, endogenous participation, exogenous entry, free entry, incumbent-quality comparative statics, and compatibility. Strongest platform-theory threat, but free entrants are an anonymous symmetric fringe and network effects are platform/group-specific rather than a shared residual incumbent outside option.
- Tan and Zhou (2021): network externalities can reverse standard competition/entry effects and free entry may be excessive or insufficient. This kills the generic `network effects reverse entry effects` claim, not the network-location necessity theorem.
- Liu, Rivadeneyra and Reshidi (2026): public entry can reallocate users away from a private network and create network fragmentation; only one public and one private platform are modeled, so no bilateral entrant cross-effect is defined.
- Economides (1996), Amir and Lazzati (2011), Gama (2019), Gama et al. (2020), Gama and Samano (2021), Basak (2021), Church and Gandal (1996), and related installed-base literature: network effects/compatibility materially alter entry, viability, and welfare, but the audited models do not establish the candidate residual-core bilateral necessity/sufficiency result.
- Common-enemy literature: stronger common threats can transform cooperation games, so `a stronger common rival induces cooperation` is not available as a novelty claim. The candidate mechanism instead works through endogenous reallocation that changes the common residual outside option.

### Stage-6 contribution narrowing

The following broad claims are killed or demoted:

- `network-effect location matters for entry` — **KILLED AS TOO BROAD** by compatible versus firm-specific network theory;
- `network effects can reverse entry effects` — **KILLED AS PRIOR ART**;
- `a stronger common rival induces coordination` — **KILLED AS PRIOR ART / TOO BROAD**;
- anti-coordination to coordination — **COROLLARY / EQUILIBRIUM IMPLICATION**;
- excessive proliferation — **WELFARE COROLLARY**;
- public/regional ownership — **APPLICATION ONLY FOR THE CENTRAL THEOREM**.

The surviving candidate is narrower:

> In a bilateral entry game against a residual incumbent network that serves as a shared endogenous outside option, entry by one peripheral platform can reduce the core's installed base and thereby lower the outside-option value faced by the other entrant. Core-side network feedback can therefore reverse the bilateral entry externality, whereas entrant-side feedback alone cannot do so in the canonical class.

### Canonical verdict

**GO**

Reason code: **GO — DISTINCT RESIDUAL-CORE NETWORK-LOCATION ENTRY THEOREM**.

This GO is deliberately narrower than `network-location matters`. No exact proposition-level prior art or less-than-one-page specialization was identified for the residual-core/shared-outside-option necessity result.

### Routing

`GO TO STAGE 7 — WELFARE / GENERALITY / INSTITUTIONAL VALIDATION`.

### Next-stage contract

Stage 7 may:

1. complete welfare analysis using the frozen model;
2. compare decentralized and social entry;
3. validate the institutional interpretation of the essential shared-core primitive;
4. test generality across nearby applications without changing the theory;
5. derive testable predictions where possible.

Stage 7 may not add players, pricing, compatibility, dynamics, local spillovers, a second core, or new network primitives. A substantive theory change requires rollback.

---

## Canonical Stage 7 — Welfare / Generality / Institutional Validation

### Provenance

- Stage 6 head SHA: `34a955faa1814ac8312051b971c1f75d59aad1ee`
- Stage 7 branch: `shared-core/stage7-welfare-generality-institution`
- Theory frozen: yes
- New primitives added: none

### Welfare reconstruction

Aggregate welfare was re-derived from individual utility. There is no separate consumer population, so no separate consumer-surplus object is created.

Verified net welfare:

- `SW00 = 2 C0`
- `SW10 = 2 C1 + [h^2 + (h-tau)^2]/2 - F`
- `SW11 = 2 C2 + s^2 - 2F`

Define social gross marginal benefits `G0` and `G1`. The cross-region welfare effects satisfy

`G0 = D0 + X0`

and

`G1 = D1 + X1`.

The second-entry effect on the already-entered region is

`X1 = - beta tau [2 - 2(A_L-A_M) + beta tau] / [2(1-2 beta)^2] < 0`

throughout the canonical interior domain.

### Central welfare identity

Definitions imply

`X1 - X0 = Gamma`,

hence

`G1 - G0 = 2 Gamma`.

Therefore the same `A*` that reverses private bilateral entry interaction also changes social marginal entry benefits from decreasing to increasing.

### Welfare results

1. **Efficient anti-coordination exists on a nonempty open set.** Strict examples exist where the asymmetric one-entry Nash equilibria coincide with `N_SP=1`.
2. **Coordination can generate socially excessive proliferation.** In every canonical coordination-multiplicity region, `(1,1)` is socially dominated by `(0,0)`.
3. **No full under-entry coordination failure exists.** There is no region in which no entry is a Nash equilibrium while `N_SP=2`.
4. **First-entry under-provision can occur.** In the substitutes region, if the first-entry cross-region effect `X0` is positive, there is an interval with no-entry Nash equilibrium but one-entry social optimum.

### Generality

The required residual-core structure was tested across multiple environments.

Strong plausible non-platform mappings were found in:

- peripheral shopping centers versus a common downtown retail agglomeration, where shopper diversion can thin downtown footfall/variety agglomeration;
- regional airports versus a dominant hub, where hub-bypass traffic can thin a common hub whose connectivity value depends on concentrated traffic and connections.

Labor-market channels are plausible but not exact because multi-channel search is common. Local marketplaces are near-exact but remain within the platform family. Technology-transfer and co-creation intermediaries more naturally exhibit complementarity to national/metropolitan ecosystems.

### Institutional validation of co-creation hubs

- Shared metropolitan/national core: **SUPPORTED**.
- Core value from information/human-network thickness: **SUPPORTED**.
- Common core across multiple regions: **SUPPORTED**.
- Cross-regional local-hub friction: **PLAUSIBLE BUT INDIRECT**.
- Binary establishment as reduced form: **SUPPORTED**.
- Regional hub entry materially diverts participation from the shared core: **PLAUSIBLE BUT INDIRECT / WEAK**.

The last item is the binding blocker. Ehime's official 2026 metropolitan co-creation project aims to strengthen E:N BASE members' network with metropolitan firms and create mutual collaboration. Ehime's DX plan treats digital, metropolitan and local hubs as a portfolio. NEXs Tokyo likewise helps regional startups use metropolitan resources. This does not prove zero substitution in participant time or attention, but it prevents treating residual-core depletion as an established institutional fact for co-creation hubs.

### Co-creation-hub application verdict

**WEAK FIT**.

The general theory survives. The original application is not yet strong enough to support a literal institutional interpretation of the theorem.

### Predictions

Externally meaningful predictions survive concerning:

- the sign of neighbor-entry interaction as core attractiveness changes;
- joint-entry clustering across the threshold;
- stronger reversal where core value is more installed-base-sensitive;
- lower strategic threshold when cross-peripheral friction rises;
- the identical private/social interaction threshold implied by `G1-G0=2 Gamma`.

### Canonical Stage-7 verdict

**CONDITIONAL GO**

Reason code: **CONDITIONAL GO — REGIONAL-TO-CORE DIVERSION EVIDENCE WEAK**.

Welfare gate: PASS.

Generality gate: PASS.

Prediction gate: PASS.

Institutional gate for the motivating application: one unresolved blocker.

### Authorized next action

Resolve exactly one question:

> Does participation in regional co-creation/innovation hubs substitute at a meaningful margin for direct participation in the shared metropolitan/national ecosystem, such that regional establishment can thin the residual core faced by firms in other regions?

No theory modification is authorized.

Do not proceed to Stage 7.5 while co-creation hubs remain the intended headline application and this blocker remains unresolved. If the blocker cannot be supported, formally demote/replace the application before Full-Theory Freeze rather than altering the model.