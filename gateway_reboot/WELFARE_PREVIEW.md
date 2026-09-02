# Welfare Preview — Stage 4R

This is a preview only; no policy section is authorized.

Gross regional utilities are `r^2/2`, `h^2/2`, `x^2/2`, `s^2/2` across the four entry states. Aggregate net welfare is

- `SW00=r^2`;
- `SW10=(h^2+x^2)/2-F`;
- `SW11=s^2-2F`.

Define social gross marginal benefits

`G0=SW10(F=0)-SW00`,

`G1=SW11(F=0)-SW10(F=0)`.

Regional entry incentives satisfy

`G0=D0+X0`, `G1=D1+X1`,

where the cross-region effects are

`X0=(x^2-r^2)/2 = [rho-(1-beta)tau][2B+rho-(1-beta)tau]/[2q^2] >0`,

`X1=(s^2-h^2)/2 = beta tau[2B+2rho-beta tau]/[2q^2] >0`.

Thus each regional government under-internalizes a positive effect on the other region.

Exact identities:

`X1-X0=Gamma^G`,

`G1-G0=2 Gamma^G`,

`(G0+G1)/2 = D0+X1 = D1+X0 > max{D0,D1}`.

## Coordination region

When `Gamma^G>0`, `D0<D1`. For `D0<F<D1`, both `(0,0)` and `(1,1)` are Nash equilibria. Yet full entry is the unique social optimum:

- `SW11-SW10=G1-F>D1-F>0`;
- `SW11-SW00=G0+G1-2F=2[(G0+G1)/2-F]>2(D1-F)>0`.

Hence the coordination failure is **under-provision**, not excessive duplication: no entry is a decentralized equilibrium although both gateways are socially optimal.

## Substitutes region

When `Gamma^G<0`, `D1<D0`. For `D1<F<D0`, asymmetric one-gateway profiles are Nash equilibria. Since `G0>D0>F`, zero entry is socially inferior. The planner chooses one gateway if `F>G1` and two if `F<G1`; decentralized anti-coordination can therefore be efficient or still under-provide.

No over-provision result appears in the canonical model because both cross-region effects are strictly positive.