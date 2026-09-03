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

The slight difference between the finite-difference production stationary point
and the exact algebraic public-FOC root is expected because the production
solver identifies the public FOC by finite differences. The exact eliminant is
the authority for this L3-1 algebraic branch.

## Branch structure

The degree-105 projection factors into multiple components. Hence root ordering
and component selection are substantive parts of the theorem; a single
unqualified resultant equation would include extraneous or non-target branches.

The next required geometry step for a full Level-3 theorem is to pair the G
algebraic price component with the distinct B3 equilibrium and then impose all
regular-interior and SOC inequalities on the paired root.
