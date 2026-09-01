# Stage 4 — Proposition Audit

Date: 2026-09-01

| Proposition | Algebraically true? | Requires endogenous participation? | Nonempty region? | Prior-art distinct? | Verdict |
|---|---:|---:|---:|---:|---|
| P1: `alpha_R > alpha_SP` and `n_R < n_SP` | Yes on `Omega` | Partly. The exact interior ordering uses participation, but a planner/regional wedge can persist under fixed participation through outsider welfare weights | Yes | No | **KILL AS CONTRIBUTION** |
| P2: interior participation threshold/tipping | No | N/A | No | N/A | **KILL** |
| P3: excessive home bias eventually lowers the regional sponsor's own payoff | Yes | **Yes**. With fixed `n`, regional payoff is strictly increasing in `alpha` | Yes | No convincing theorem-level distance | **MATHEMATICALLY SURVIVES; NOVELTY FAILS** |
| P4: neutrality commitment improves sponsor payoff | Not separately developed | Would depend on P3 | N/A | Weak / standard commitment logic | **DO NOT PURSUE** |

## P1 — Composition wedge

On

`Omega={0<v<2, B>1/2, v/2+v(B-1/2)/sqrt(1+4v)<A<vB}`,

both objectives are strictly concave and interior.

At `alpha_R`,

`U_R'(alpha_R)=0`,

while

`W'(alpha_R)=n(alpha_R)n'(alpha_R)<0`.

Hence

`alpha_SP<alpha_R`

and

`n_SP>n_R`.

However, with participation fixed, a planner/regional difference can still arise because the planner directly counts outsider access utility while the regional sponsor does not. Therefore P1 is not a clean mechanism-essential contribution.

**Verdict: KILL AS CONTRIBUTION.**

## P2 — Participation threshold/tipping

The equilibrium participation function is

`n(alpha)=[sqrt(1+4v(1-alpha))-1]/2`.

It is positive and smooth for every `alpha<1` and reaches zero only at `alpha=1`.

Creating an interior threshold would require an additional fixed cost, minimum cost, discrete participation mass, coordination effect, or similar extra primitive. These are not authorized and would violate the no-feature-accumulation rule.

**Verdict: KILL.**

## P3 — Local-welfare reversal

Regional payoff is

`U_R(alpha)=A alpha+B n(alpha)`.

On the regional interior region,

`U_R''(alpha)<0`,

`U_R'(0)>0`,

and

`U_R'(1)<0`.

Therefore a unique `alpha_R in (0,1)` maximizes regional own welfare, and

`U_R'(alpha)<0` for all `alpha>alpha_R`.

This is mechanism-essential mathematically. With fixed participation,

`U_R^F(alpha)=A alpha+B nbar`

and

`dU_R^F/dalpha=A>0`,

so the reversal disappears.

Nevertheless, after label stripping the proposition is:

> A biased principal that values both favoritism and the size of an endogenously participating user pool optimally limits its bias because excessive bias drives users away.

That is too close in economic structure to established platform/intermediary-bias participation trade-offs. Public regional sponsorship does not add a new strategic term.

**Verdict: MATHEMATICALLY SURVIVES; NOVELTY FAILS.**

## P4 — Neutrality commitment

The regional sponsor already internalizes the effect of bias on its own participant-pool value and chooses `alpha_R` optimally. A stricter exogenous cap below `alpha_R` would lower regional payoff; a nonbinding cap adds nothing.

Obtaining a distinct commitment theorem would require a separate commitment failure, agency problem, time inconsistency, or governance game. Those are new mechanisms and prohibited in Stage 4.

**Verdict: DO NOT PURSUE.**

## Overall proposition verdict

The minimal model passes the mathematical gate only for P1 and P3. P2 fails and P4 has no independent content.

P3 is the strongest internal result, but it does not survive the novelty/stripping gate as a defensible contribution.

**Stage-level implication: NO-GO — RESULT TOO WEAK.**