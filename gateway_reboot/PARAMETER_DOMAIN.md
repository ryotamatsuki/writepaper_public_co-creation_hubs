# Canonical Parameter Domain and Strategic Reversal

Let `q=1-2 beta` and `d=rho-tau`.

## Canonical interior domain

A clean domain supporting the intended gateway structure is

`0<beta<1/2`,

`0<tau<rho<q`,

`0<B<q-rho`.

Then `0<r<x<h<s<1`, cross-region gateway use is strictly better than direct access when only the neighbor gateway exists, and the fixed point is unique.

## Threshold

Because `d Gamma^G/dB = -d/q^2<0`, there can be at most one crossing in `B`.

The exact root is

`B* = beta(1-beta) tau^2/d - d/2`.

Thus, whenever `B*` lies inside `(0,q-rho)`,

- `B<B* => Gamma^G>0` (strategic complements),
- `B>B* => Gamma^G<0` (strategic substitutes).

An exact convenient condition for an interior crossing is

`d^2 < 2 beta(1-beta) tau^2 < d[2(q-rho)+d]`.

Equivalently, the first inequality gives `B*>0` and the second gives `B*<q-rho`. These are strict inequalities and therefore define a nonempty open set.

## Open-set calibration

Take

`beta=0.10`, `tau=0.20`, `rho=0.23`.

Then `q=0.80`, `d=0.03`, the canonical `B` domain is `(0,0.57)`, and

`B*=0.105`.

At `B=0.05`:

`(r,x,h,s)=(0.0625,0.125,0.325,0.35)`,

`D0=0.050859375`, `D1=0.0534375`, `Gamma^G=+0.002578125`.

At `B=0.20`:

`(r,x,h,s)=(0.25,0.3125,0.5125,0.5375)`,

`D0=0.100078125`, `D1=0.095625`, `Gamma^G=-0.004453125`.

All inequalities are strict, so continuity gives open neighborhoods with opposite signs.

## Comparative statics of the threshold

Holding the other primitives fixed:

`dB*/dbeta = (1-2beta)tau^2/(rho-tau) > 0`.

`dB*/drho = -beta(1-beta)tau^2/(rho-tau)^2 - 1/2 < 0`.

`dB*/dtau = beta(1-beta)[2tau/(rho-tau)+tau^2/(rho-tau)^2] + 1/2 > 0`.

Interpretation: stronger shared-core network feedback or greater cross-region friction makes complementarity persist up to a stronger baseline core; stronger own-gateway access improvement relative to cross-region friction increases overlap and moves the system toward substitution earlier.