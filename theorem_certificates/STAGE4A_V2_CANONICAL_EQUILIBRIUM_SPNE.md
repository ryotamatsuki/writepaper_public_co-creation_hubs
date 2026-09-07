# Theorem Certificate — Canonical Public Equilibrium / SPNE Status

Claim ID: canonical G/B3 equilibrium status underlying `prop:exactcertificate`, model Section 2, and welfare interpretation.

Stage-4A date: 2026-09-07

## Claim being audited

The manuscript states that the game is studied in subgame-perfect equilibrium and refers to the canonical G/B3 algebraic stationary pair as public/full-game equilibria. The primitive public strategy set is `x_i in [0,1]`; the private strategy set is `p_T>=0`; downstream project route choice includes `0,H1,H2,HT` and partner participation follows the stated primitive distribution.

For the G canonical pair to be a public Nash/SPNE outcome, each government must have no profitable finite deviation after a valid downstream private-price and participation continuation is re-solved. B3 likewise requires no profitable finite public deviation at the fixed matched price if it is called a public Nash equilibrium.

## Production evidence

The exact Krawczyk/interval certificate establishes a unique coupled stationary root inside a narrow rational isolating box, strict local curvature, nonsingularity, and signs of local cross derivatives. It does not certify global public best responses over `x_i in [0,1]` or all downstream participation regimes.

The old Stage-4 numerical routine explicitly solves only the `0 -> H_T -> H_i` regime and assigns `-1e9` private profit when that regime solver returns invalid. Under v2.0, this cannot be used as a global continuation or deviation certificate.

## Independent reconstruction

FAIL for global equilibrium; PASS for existence of the audited off-regime continuation used to falsify it.

A clean-room evaluator reconstructs project choice from the primitive upper envelope, clips partner participation according to `c~U[0,1]`, solves the participation fixed point from multiple starts, and re-optimizes the private price on the full positive-profit interval.

## Finite-deviation result

At canonical G, government 1 can deviate from `x_1 ~= 0.68402819636` to `x_1=0.1875` with `x_2` fixed. After downstream re-optimization:

- candidate direct `W_1 ~= 0.22301946658`;
- deviation direct `W_1 ~= 0.25444926939`;
- gain `~= +0.03142980281`.

The deviation activates rival-public use and makes the deviator's own public hub inactive. It is therefore outside the production central regime but inside the stated primitive strategy/allocation domain.

At canonical B3 the same finite deviation raises regional welfare by approximately `+0.02795444338`.

## Continuation completeness

FAIL.

The previous workflow did not provide valid continuations for all material off-path public deviations and regime switches. The counterexample demonstrates that this omission is economically material rather than merely formal.

## Boundary/corner/regime audit

FAIL.

A profitable rival-public-entry / own-public-inactivity regime was omitted from the global best-response check.

## Benchmark-definition status

B3 remains a valid matched-price diagnostic as a local stationary benchmark. It may not be called a public Nash equilibrium at the canonical vector without a repaired global solution.

## Manuscript-scope consequence

The following wording is not certified:

- global `subgame-perfect equilibrium` for the canonical G pair;
- `public equilibrium` / `full equilibrium` for the canonical stationary pair;
- welfare conclusions stated as properties of the decentralized equilibrium when evaluated only at that stationary pair.

The exact interval certificate itself remains valid as a sign-at-stationary-algebraic-root certificate.

## Certificate state

`FAIL — PROFITABLE FINITE OFF-REGIME DEVIATION`.

Earliest affected stage: **Stage 4**.

Per v2.0 routing, Stage 8 theory freeze and all dependent Stage 9--14 outputs are stale until the game is globally re-solved or the equilibrium/contribution object is formally redefined through the appropriate earlier-stage change control and then recertified.