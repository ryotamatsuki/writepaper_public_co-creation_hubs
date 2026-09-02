# Gateway Participation Equilibrium

Define `q=1-2 beta` and assume the canonical interior domain `0<beta<1/2`, `0<tau<rho<q`, `0<B<q-rho`.

For a region with best available gateway bonus `a_i`, the connected cutoff is

`z_i = B + beta N + a_i`,

with `N=z_1+z_2`. Hence the common fixed point is unique because the interior aggregate slope is `2 beta<1`.

## No gateways `(0,0)`

Both regions use direct access, `a_1=a_2=0`:

`r = B/q`, `N_00=2r=2B/q`.

## One gateway `(1,0)`

Region 1 uses its own gateway, bonus `rho`. Region 2 uses region 1's gateway, bonus `rho-tau`; direct access is strictly dominated because `rho>tau`.

`h = (B+rho-beta tau)/q`,

`x = h-tau = [B+rho-(1-beta)tau]/q`,

`N_10=h+x=(2B+2rho-tau)/q`.

## Two gateways `(1,1)`

Each region uses its own gateway; own access dominates the other-region gateway by `tau`:

`s=(B+rho)/q`, `N_11=2s=2(B+rho)/q`.

## Ordering and feasibility

Under the canonical domain,

`0<r<x<h<s<1`.

Useful differences:

- `x-r=[rho-(1-beta)tau]/q>0`;
- `h-x=tau>0`;
- `s-h=beta tau/q>0`.

Thus it is sufficient and, within `0<tau<rho`, equivalent for the intended interior configuration to require `B>0` and `B+rho<q`.

Gateway entry always thickens the shared core:

- `N_10-N_00=(2rho-tau)/q>0`;
- `N_11-N_10=tau/q>0`.

The second increment raises every connected user's shared-core value by `beta tau/q`.