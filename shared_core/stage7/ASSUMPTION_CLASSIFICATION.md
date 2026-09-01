# Stage 7 — Essential vs Tractability Assumptions

| Assumption | Essential to mechanism? | Tractability only? | Institutional/theoretical defense | Stage-7 status |
|---|---|---|---|---|
| Two peripheral entrants | YES | No | bilateral cross-entry externality is the object of interest | ESSENTIAL |
| Shared residual core | YES | No | without a common outside option there is no residual-core feedback | ESSENTIAL |
| Core-side network feedback | YES in canonical class | No | Proposition 2: entrant-side effects alone do not generate complementarity; network-valued cores have clear platform/retail/airport analogues | ESSENTIAL |
| Endogenous user allocation | YES | No | entry must change core participation for the mechanism to operate | ESSENTIAL |
| Single-homing | Essential to the current closed-form canonical proof; not established as generally essential | Partly | defensible per transaction/episode in retail or travel; weak as a literal description of startup ecosystems or labor search | MATERIAL LIMITATION |
| Uniform \(c\) | No | YES | Stage 4 authorized robustness with \(F(c)=c^\eta\) near \(\eta=1\) preserves the reversal | TRACTABILITY |
| Linear core network effect | No | YES | Stage 4 authorized concave-core robustness preserves the reversal locally | TRACTABILITY |
| Regional symmetry \(A_{L1}=A_{L2}\) | No | YES | Stage 4 small-asymmetry check preserves both regions' reversal | TRACTABILITY |
| Binary entry | Essential to the exact anti-coordination/coordination game-form result; probably not to the broad feedback intuition | Partly | large fixed establishment decisions are naturally discrete in facilities/infrastructure | ESSENTIAL TO CURRENT GAME FORM / REDUCED FORM |
| \(\tau>0\) | YES for distinct peripheral options and cross-use friction | No | geographic/relationship/access frictions are economically meaningful; exact magnitude is application-specific | ESSENTIAL |
| Public/regional welfare objective | NO for the central strategic reversal | No | Stage 3–4 diagnostics show a simple entrant objective can exhibit the residual-core reversal | NONESSENTIAL / APPLICATION |
| \(\beta_L=0\) baseline | NO | YES/minimality | core-only feedback is sufficient; local-only feedback is insufficient; positive local feedback belongs only in the diagnostic robustness model | MINIMALITY CHOICE |
| Fixed establishment cost \(F\) | YES for entry regime selection | No | large discrete facilities/programs/market entry require resource commitments | ESSENTIAL TO ENTRY GAME |
| Identical shared core for both regions | YES in canonical class | No | national platforms, downtown cores or hub airports can provide a materially overlapping common option | ESSENTIAL |

## Essential causal chain

The mechanism can be stated without canonical notation:

1. two peripheral options consider entry;
2. both compete against the same network-dependent incumbent/core;
3. entry by one peripheral option diverts users from the core;
4. the thinner core becomes less valuable to users in both peripheral markets;
5. this weakens the outside option faced by the other entrant;
6. when the residual-core effect dominates direct business stealing, bilateral entry becomes complementary.

Steps 2–4 are the truly essential primitives. Uniform heterogeneity, linearity, symmetry and zero local network feedback are simplifying devices.

## Existing authorized alternative formulations

No new Stage-7 robustness model was introduced. The following prior authorized checks are carried forward:

- **Concave core network:** reversal survives for nearby \(u_M=A_M+\beta(n_M-\kappa n_M^2)\).
- **Small regional asymmetry:** reversal survives when \(A_{L1}=A_L+\epsilon\), \(A_{L2}=A_L-\epsilon\) for small \(\epsilon\).
- **Alternative cost distribution:** reversal survives for \(F(c)=c^\eta\) with \(\eta\) near 1.

These are sufficient for Stage 7 to classify linearity, symmetry and uniformity as tractability assumptions rather than essential mechanisms.

## Main unresolved scope limitation

Single-homing remains the largest untested theoretical restriction. The current workflow prohibits adding multi-homing at Stage 7. Therefore the project must not claim that the theorem is robust to multi-homing. This is a scope limitation to be considered explicitly at Stage 7.5, not an authorization for theory expansion.
