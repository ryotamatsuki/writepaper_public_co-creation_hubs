# Stage 4 — Model Definitions

Date: 2026-09-01

Stage 4 base SHA: `a4bff663b34898fd61a7eb67d4ac989515b5aca5`

Workflow SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`

Branch: `research-reboot/stage4-minimal-neutrality-model`

## 1. Purpose

This file defines the smallest model permitted by the Stage 3 contract. It is a kill-test representation of M10, not a proposed final paper model.

The only substantive mechanism is:

`home bias -> lower expected outsider access -> lower outsider participation -> thinner outside option set -> lower regional value from intermediation`.

No legacy hub/pipeline primitive is used.

## 2. Players and timing

1. A regional sponsor `R` chooses a home-priority parameter `alpha in [0,1]`.
2. A unit mass of potential outside participants observes `alpha` and decides whether to enter the intermediary's searchable pool.
3. The intermediary implements the priority rule under fixed matching/attention capacity. No independent intermediary optimization is introduced.
4. A social planner is a benchmark only and chooses the same `alpha` under the same participation technology.

## 3. Outside participation

Potential outside participants have participation costs

`c ~ Uniform[0,1]`.

Let `n in [0,1]` be the endogenous mass of outside participants.

The intermediary has fixed outside-facing attention. Given home bias `alpha`, an outsider's gross expected access/match value is

`g(alpha,n) = v(1-alpha)/(1+n)`,

where `v>0` is the gross value of access and `1/(1+n)` is congestion from a thicker applicant pool under fixed attention capacity.

An outsider participates iff

`c <= g(alpha,n)`.

For an interior cutoff under the uniform distribution,

`n = g(alpha,n)`,

so equilibrium participation solves

`n(1+n) = v(1-alpha)`.

Hence

`n(alpha) = [sqrt(1+4v(1-alpha))-1]/2`.

To keep the uniform-CDF cap inactive for every `alpha`, use the benchmark restriction

`0 < v < 2`.

Then `n(0)<1`, `n(alpha)>0` for every `alpha<1`, and `n(1)=0`.

## 4. Economic interpretation of the regional payoff

The sponsor obtains two forms of local benefit.

First, greater home priority directly raises local-favored intermediation value at constant marginal value `A>0`.

Second, the sponsor obtains local option/search value from having a thicker qualified outside pool, at marginal value `B>0`.

The reduced-form regional payoff is therefore

`U_R(alpha,n) = A alpha + B n`.

The term `B n` should be read as the local value of breadth/quality in the outside opportunity set assembled by the intermediary, not as outsider welfare. The model deliberately does not insert a separate cross-region spillover coefficient.

This specification is useful for the Stage 4 kill test because any local-welfare reversal must come from endogenous pool thinning. With participation artificially fixed, the regional payoff is strictly increasing in `alpha`.

## 5. Outsider surplus and planner welfare

At an endogenous participation equilibrium the marginal entrant has cost `c*=n`, because the participation condition implies gross entrant value `g=n`.

Aggregate outside net surplus is therefore

`S_X = integral_0^n (n-c) dc = n^2/2`.

The social planner uses exactly the same participation and priority technology and maximizes

`W(alpha,n) = A alpha + B n + n^2/2`.

The planner is not given a better technology.

## 6. Backward-induction objects

Define

`D(alpha) = 1 + 4v(1-alpha)`.

Then

`n(alpha) = [sqrt(D)-1]/2`,

`dn/dalpha = -v/sqrt(D) < 0`,

and

`d2n/dalpha2 = -2v^2/D^(3/2) < 0`.

The reduced regional objective is

`U~_R(alpha) = A alpha + B n(alpha)`.

The reduced planner objective is

`W~(alpha) = A alpha + B n(alpha) + n(alpha)^2/2`.

## 7. Candidate interior parameter region

Let

`s = sqrt(1+4v)`.

The open region used for the joint interior solution is

`Omega = { 0<v<2, B>1/2, v/2 + v(B-1/2)/s < A < vB }`.

On `Omega`:

- regional and planner objectives each have a unique interior optimum;
- outside participation is interior and below one;
- the regional sponsor chooses more home bias than the planner;
- excessive home bias beyond the regional optimum reduces the sponsor's own payoff.

The derivations and kill tests are recorded separately.

## 8. Fixed-participation counterfactual

To test whether participation is mechanism-essential, fix outside participation at `nbar in (0,1)` while retaining the same priority/access technology.

Regional payoff becomes

`U_R^F(alpha) = A alpha + B nbar`,

so

`dU_R^F/dalpha = A > 0`.

Thus the regional sponsor selects `alpha=1` and there is no local-welfare reversal.

For welfare accounting, if the same mass `nbar` remains in the pool, outsider gross access value is

`v(1-alpha)/(1+nbar)`

and outsider net surplus is

`nbar*v(1-alpha)/(1+nbar) - nbar^2/2`.

This fixed-participation benchmark is used to separate a genuine participation feedback from a simple welfare-weight difference.

## 9. What this model does not contain

No second strategic region, intermediary competition, private intermediary, dynamics, absorptive capacity, KPI contracting, political capture, endogenous network formation, fixed hub cost, R&D stage, migration, or legacy hub/pipeline technology is present.