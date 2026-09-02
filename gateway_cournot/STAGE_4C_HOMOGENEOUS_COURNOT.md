# Stage 4C — Homogeneous Cournot Gateway Hard Kill

## Executive verdict

**GO — HOMOGENEOUS COURNOT PRESERVES ENDOGENOUS GATEWAY STRATEGIC REVERSAL.**

Product differentiation is **NOT NECESSARY**.

## Complete game

1. Governments choose `E_i in {0,1}`.
2. Firms choose ecosystem engagement `x_i>=0`.
3. Firms choose homogeneous Cournot quantities `q_i>=0` under `P=a-q_1-q_2`.

Gateway provision lowers the firm's quadratic engagement-cost coefficient. Engagement lowers marginal cost directly and through a common network term:

`c_i=cbar-x_i-beta(x_i+x_j)`.

## Cournot equilibrium

Let `A=a-cbar`, `m=2+beta`, `n=beta-1`.

For given engagement,

`q_i^*=[A+m x_i+n x_j]/3`.

Thus

`dq_i^*/dx_j=(beta-1)/3`.

Operating profit is `(q_i^*)^2`.

## Engagement equilibrium

For profile-specific engagement costs `(K_i,K_j)`, define

`a_i=9K_i-2m^2`, `a_j=9K_j-2m^2`, `b=2mn`, `Delta=a_i a_j-b^2`.

Then

`x_i^*=2mA(a_j+b)/Delta`,

`q_i^*=3A K_i(a_j+b)/Delta`,

`pi_i^*=A^2 K_i a_i(a_j+b)^2/Delta^2`.

The engagement best-response slope is `b/a_i`, so engagement is substitutable for `beta<1` and complementary for `beta>1` under the stable domain.

## Stability domain

A clean sufficient domain is

`A>0`, `beta>=0`, `0<tau<rho<kappa`, `F>0`,

and

`kappa-rho > [2(2+beta)^2 + |2(2+beta)(beta-1)|]/9`.

This gives strict own concavity, unique joint engagement equilibrium and positive engagement/output for every gateway profile.

## Government objective

Canonical baseline:

`W_i=pi_i+CS/2-FE_i`.

Producer-only `pi_i-FE_i` is retained as an instructed comparison.

Define `H=kappa`, `L=kappa-rho`, `M=kappa-rho+tau`, and let `R_i(K_i,K_j)=pi_i+CS/2` after solving the complete firm subgame.

`D_0^C=R_i(L,M)-R_i(H,H)`,

`D_1^C=R_i(L,L)-R_i(M,L)`,

`Gamma^C=D_1^C-D_0^C`.

The exact rational function is high-order and does not admit a useful global economic factorization. Mechanism identification is instead established through exact deletion results and strict open-set witnesses.

## Network necessity

For the regional-welfare baseline, `beta=0` implies

**`Gamma^C<0` everywhere in the canonical stable domain.**

The proof rewrites the exact expression using `ell=kappa-rho`, `d=rho-tau`, `t=tau`, `y=ell-4/3`; the numerator after extracting the leading minus sign is a polynomial with strictly positive coefficients in `d,t,y`.

Thus common-network feedback is necessary for gateway complementarity under the canonical public/regional objective.

## Both signs on open sets

With `A=1, kappa=5, rho=3/2, tau=6/5`:

- at `beta=1/2`, `Gamma^C=-27845669343/2525102237378<0`;
- at `beta=6/5`, `Gamma^C=7916917186023346560/126835784788527792481>0`.

All stability inequalities are strict. Continuity gives nonempty open sets. Numerical perturbation checks preserve the negative sign under +/-5% perturbations of every primitive and the positive sign under +/-2% perturbations.

## A is not a strategic-sign threshold

All gross equilibrium payoffs scale with `A^2`. Therefore `A=a-cbar` changes the level of gateway willingness to pay but not `sign(Gamma^C)`. The old Stage 4R baseline-core threshold does not survive explicit Cournot competition.

## Government-objective comparison

Producer-focused government also has both signs, so profit objective alone is sufficient for mathematical reversal. However, at `beta=0` it can already yield positive `Gamma^C` due pure rent shifting/business stealing. Regional welfare changes the interpretation qualitatively and gives the clean network-necessity theorem.

Classification: **RESULT CHANGES QUALITATIVELY**.

## Welfare preview

National welfare is `pi_1+pi_2+CS-F(E_1+E_2)`. Product-market rivalry makes cross-region externalities sign-ambiguous. Both decentralized under-provision and over-provision occur in admissible parameter regions.

## Numerical hard kill

500,000 admissible draws were generated with the sufficient stability condition built into sampling.

Producer objective:

- `Gamma^C>0`: 42,831
- `Gamma^C<0`: 457,169
- near zero: 0

Regional-welfare objective:

- `Gamma^C>0`: 41,371
- `Gamma^C<0`: 458,629
- near zero: 0

Independent direct solutions of the engagement linear system were compared with the closed forms for 10,000 draws across all four gateway profiles:

- maximum engagement error: `2.67e-15`;
- maximum engagement-FOC residual: `1.07e-14`.

No product-differentiation parameter was used.

## Institutional verdict

**PLAUSIBLE FIT.**

The model represents co-creation hubs as access infrastructure linking local firms to a wider ecosystem, while firms remain genuine product-market competitors. This matches the direction of current E:N BASE and metropolitan-linkage institutional evidence better than the earlier core-depletion model.

## Previous model classification

**PARTIAL DIAGNOSTIC ONLY.**

## Final routing

**FREEZE HOMOGENEOUS-COURNOT GATEWAY MODEL.**

Then proceed to an exact proposition-level prior-art re-kill. Do not introduce differentiated Cournot before that re-kill.