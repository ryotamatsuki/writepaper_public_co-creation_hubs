# Stage 4R — Minimal Gateway Model Falsification

## Executive verdict

**GO — INTERIOR GATEWAY STRATEGIC REVERSAL SURVIVES MINIMAL MODEL.**

No unauthorized primitive was required.

## Final minimal model

Let `B=A_M-delta`; retain only `(B,beta,rho,tau,F)` with firm cost `c~U[0,1]` and outside option 0.

Core value is `B+beta N`. Access-route bonuses are 0 (direct), `rho` (own gateway), and `rho-tau` (other gateway). Regional government maximizes origin-region realized firm utility net of `F E_i`.

Canonical interior domain:

`0<beta<1/2`, `0<tau<rho<1-2beta`, `0<B<1-2beta-rho`.

## Participation cutoffs

With `q=1-2beta`:

`r=B/q`,

`h=(B+rho-beta tau)/q`,

`x=h-tau=[B+rho-(1-beta)tau]/q`,

`s=(B+rho)/q`.

They satisfy `0<r<x<h<s<1`.

## Entry incentives

`D0=(rho-beta tau)(2B+rho-beta tau)/(2q^2)`.

`D1=tau(1-beta)(2B+beta tau+2rho-tau)/(2q^2)`.

Therefore

`Gamma^G = [2 beta(1-beta) tau^2 - (rho-tau)(2B+rho-tau)]/(2q^2)`.

## Proposition — Gateway strategic reversal

Let `d=rho-tau>0`. If

`d^2 < 2 beta(1-beta) tau^2 < d[2(q-rho)+d]`,

then

`B*=beta(1-beta)tau^2/d-d/2`

lies strictly inside `(0,q-rho)`. Because

`partial Gamma^G/partial B = -d/q^2<0`,

there is a unique crossing:

- `B<B*`: `Gamma^G>0`, strategic complements;
- `B>B*`: `Gamma^G<0`, strategic substitutes.

The inequalities are strict, so the result holds on a nonempty open parameter set.

## Mechanism identification

`beta=0` gives `Gamma^G=-d(2B+d)/2<0`: shared-core network feedback is necessary for complementarity in this canonical class.

At `tau=rho` (limit from the overlap domain), other-region gateway access offers no incremental direct-access advantage and the overlap term vanishes; `Gamma^G=beta(1-beta)rho^2/q^2>0`.

At `tau=0`, gateways are perfectly interchangeable and `Gamma^G=-rho(2B+rho)/(2q^2)<0`.

Thus shared-core thickening and gateway overlap exert opposite strategic effects.

## Comparative statics

`partial B*/partial beta=(1-2beta)tau^2/(rho-tau)>0`.

`partial B*/partial rho=-beta(1-beta)tau^2/(rho-tau)^2-1/2<0`.

`partial B*/partial tau=beta(1-beta)[2tau/(rho-tau)+tau^2/(rho-tau)^2]+1/2>0`.

## Public-objective relevance

Natural non-price private/user-count objectives yield fixed-sign strategic effects and do not reproduce the reversal. The regional welfare objective is therefore **ESSENTIAL FOR THE SIGN-SWITCH THEOREM** in this minimal formulation.

## Welfare preview

Both interregional welfare effects are positive: `X0>0`, `X1>0`, with `G1-G0=2 Gamma^G`. In the coordination interval `D0<F<D1`, no entry and full entry are both Nash equilibria, but full entry is the unique social optimum. The gateway version therefore naturally generates decentralized under-provision.

## Numerical falsification

300,000 random draws satisfying the canonical interior domain produced:

- 54,261 cases with `Gamma^G>0`;
- 245,739 cases with `Gamma^G<0`;
- 67,448 parameter triples with `B*` strictly inside the feasible `B` interval;
- zero threshold-sign mismatches;
- `X0>0` and `X1>0` in every draw;
- `G1-G0=2 Gamma^G` to machine precision;
- no counterexample to `(G0+G1)/2>max(D0,D1)`.

A 10,000-draw independent numerical fixed-point iteration matched analytical cutoffs with maximum error below `5e-13`.

## Canonical routing

**GO TO STAGE 6R — EXACT PRIOR-ART RE-KILL OF GATEWAY ENTRY THEOREM.**

No novelty claim is authorized before Stage 6R.