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
