# Selected Platform Architecture

## Name

**Project Single-Homing × Partner Multi-Homing with Partner-Side Public Facilitation**.

## Players

Regional governments `G_1,G_2`; public platforms `H_1,H_2`; national private platform `H_T`.

## Project side F

Each region has unit mass of resident projects. Productive project value `z~U[0,1]`. Each project chooses exactly one of `0,H_1,H_2,H_T`.

Utility:

`U_rh(z)=z[v+alpha n_h^P]-kappa_rh-1{h=T}p_T`, `U_0=0`.

Local/rival public access costs satisfy `kappa_rr=kappa_L`, `kappa_rj=kappa_L+tau`; `tau` is a real coordination cost and may be zero in robustness tests.

## Partner side P

A unit mass of support/collaboration partners with real per-platform participation cost `c~U[0,1]`. Partners may multi-home. Participation on each platform is a separate activity, so a partner joins h whenever its platform-specific net benefit is nonnegative.

Public platform i partner value:

`rho + x_i + beta n_i^F - c`.

Private platform partner value:

`rho_T + beta n_T^F - c`.

Interior masses:

`n_i^P=rho+x_i+beta n_i^F`,

`n_T^P=rho_T+beta n_T^F`.

## Strategic controls

`G_i` chooses `x_i in [0,1]`; `x_i` is partner-side facilitation/recruitment effort. Cost `C_i=gamma x_i^2/2`.

`H_T` chooses one project-side access fee `p_T>=0`. Partner access is free in the baseline. Private real marginal operating cost may be normalized to zero.

## Welfare

`W_i` counts surplus of region-i resident projects plus region-i resident partner surplus across platforms minus `gamma x_i^2/2`. Private fees paid by resident projects are regional outflows.

National welfare adds both regional welfare functions and private profit, cancelling fee transfers.

## Timing

1. `G_1,G_2` simultaneously choose `x_1,x_2`.
2. `H_T` observes investments and chooses `p_T`.
3. projects and partners reach the static participation fixed point.
4. surplus is realized.

Equilibrium concept: SPNE with a selected/stable participation equilibrium where multiplicity arises; Stage 4 must characterize the relevant interior regime before invoking stability.

## Benchmarks

B1 remove H2; B2 remove HT; B3 fix pT; B4 set beta=0/fix partner response; G all endogenous.

## Stage 4 target

P2-R first; P1-R second. P3–P5 accounting only.