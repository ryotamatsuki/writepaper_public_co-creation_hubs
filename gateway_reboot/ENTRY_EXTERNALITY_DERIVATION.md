# Entry Externality Derivation

Regional gross utility for origin region `i` is the integral of realized utility over connected firms. With cutoff `z`, this is `z^2/2`.

Therefore

- `W~_i(0,0)=r^2/2`;
- `W~_i(1,0)=h^2/2`;
- `W~_i(0,1)=x^2/2`;
- `W~_i(1,1)=s^2/2`.

Define

`D0=(h^2-r^2)/2`,

`D1=(s^2-x^2)/2`.

Exact simplified forms are

`D0 = (rho-beta tau)[2B+rho-beta tau]/[2q^2]`,

`D1 = tau(1-beta)[2B+beta tau+2rho-tau]/[2q^2]`.

The bilateral cross-entry effect is

`Gamma^G=D1-D0`, hence

`Gamma^G = {2 beta(1-beta) tau^2 - (rho-tau)[2B+rho-tau]}/[2q^2]`.

This exactly matches the preliminary conjecture; no correction was required.

## Mechanism decomposition

Write `d=rho-tau>0`. Then

`Gamma^G = [2 beta(1-beta) tau^2 - d(2B+d)]/[2q^2]`.

The first numerator term vanishes at `beta=0` and is the shared-core-thickening contribution. The second vanishes when `d=0` (`tau=rho`), i.e. when the other region's gateway provides no direct access improvement and gateway overlap disappears.

Deletion tests support the economic interpretation:

- `beta=0`: `Gamma^G=-d(2B+d)/2<0` — gateways are strict substitutes;
- `tau=0`: `Gamma^G=-rho(2B+rho)/(2q^2)<0` — perfectly interchangeable gateways are substitutes;
- `tau -> rho`: `Gamma^G -> beta(1-beta)rho^2/q^2>0` — with overlap removed, shared-core thickening produces complements.

Thus the sign is not imposed; it is the difference between network thickening and overlapping access.