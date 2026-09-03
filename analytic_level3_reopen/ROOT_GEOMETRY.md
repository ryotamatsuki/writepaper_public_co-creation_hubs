# Root Geometry

## Canonical G branch

All canonical decimal primitives were converted to exact rationals.

After exact elimination, the economically relevant G branch lies on a
degree-31 univariate factor `T_31(t)`.

Sturm root counting gives exactly one real root in

`58489/100000 < t_G < 58490/100000`.

Endpoint signs are opposite.

An independent elimination onto the private price gives a degree-31 price
factor with exactly one real root in

`22746/1000000 < p_G < 22747/1000000`.

Again, endpoint signs are opposite.

These are exact rational isolation intervals, not floating-point confidence
intervals.

For cross-check only, the production numerical solver reports approximately

`t_G = 0.5848907466`, `p_G = 0.0227466816`.

The finite-difference production stationary point and the exact algebraic
public-FOC root need not coincide at all displayed digits; the exact eliminant
and the coupled interval certificate are the authority for L3-1.

## Coupled canonical branch

The later L3-1 closure step pairs the G algebraic component with the distinct B3
equilibrium at the same algebraic `p_G`. A five-polynomial rational Krawczyk
operator proves a unique coupled solution in the certified box for
`(t_G,q_G,p_G,t_B,q_B)` and therefore resolves the canonical branch-matching
problem exactly. Reconstruction preserves `x_B != x_G`.

Exact rational interval implicit differentiation on that same isolated branch
certifies the G/B3 SOC and cross-sign obligations.

## Branch structure

The degree-105 projection factors into multiple components. Hence root ordering
and component selection remain substantive parts of any non-calibrated theorem;
a single unqualified resultant equation would include extraneous or non-target
branches.

The next geometry problem for a full Level-3 theorem is no longer canonical
G/B3 pairing. It is to make the branch selection uniform over a declared
primitive parameter region and project all regular-interior, SOC, and cross-sign
conditions into primitive space, for example through CAD/regular chains or
uniform Sturm/Thom conditions.
