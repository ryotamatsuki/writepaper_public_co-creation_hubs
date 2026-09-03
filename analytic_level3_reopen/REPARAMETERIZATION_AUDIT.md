# Reparameterization Audit

## Winning coordinates: `(t,q_T,p)`

At a symmetric regular-interior state let

- `delta = alpha beta`,
- `a = kappa_T+p`,
- `d = kappa_L-a`,
- `A_T = v+alpha rho_T`.

The participation block is

`q_T^2-(A_T+2 delta t)q_T+2 delta a=0`

and the public investment is reconstructed exactly as

`x=(q_T+d/t-delta(1-t)-v)/alpha-rho`.

Thus `x` is not needed as an elimination variable.

The same equation can be inverted for the network strength:

`delta = q_T(q_T-A_T)/(2(t q_T-a))`

whenever `t q_T != a`. This is Frontier M's parameter-inversion formula.

## Jacobian coordinates

Define

`u=d/t-delta t`,

`w=1-2 delta a/q_T^2`,

`H=u w-2 delta t`.

At symmetry the participation Jacobian is

`[[u,0,-t],[0,u,-t],[-delta,-delta,w]]`

with determinant `u H`. Symmetric and antisymmetric response modes therefore
separate before any expansion.

## Other coordinates tried

- `(x,t,q_T,p)`: redundant because `x` reconstructs rationally.
- `(t,s)` with `s=a/q_T`: simplifies some admissibility inequalities but makes
  the public cost term and network inversion less transparent.
- `(u=t_1+t_2,v=t_1-t_2)`: useful for local Jacobian mode separation, but the
  private and public objectives are simplest after imposing symmetry.
- participation masses `(n_i,n_T)`: equivalent affine transformations of cutoffs;
  no lower polynomial degree was obtained.

Verdict: `(t,q_T,p)` is the preferred symmetric elimination chart.
