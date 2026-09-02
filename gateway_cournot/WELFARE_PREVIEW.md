# Welfare Preview

National gross welfare is

`SW = pi_1 + pi_2 + CS`.

For gateway profiles, define

`S_0 = SW^{10}-SW^{00}`,

`S_1 = SW^{11}-SW^{01}`.

The canonical regional government instead uses

`D_0^C = R_i^{10}-R_i^{00}`,

`D_1^C = R_i^{11}-R_i^{01}`,

with `R_i=pi_i+CS/2`.

Unlike the Stage 4R participation model, the regional-social wedges are not globally one-signed. Product-market competition can make a regional gateway confer positive or negative effects on the other region depending on network strength and rivalry.

In the 500,000-draw admissible hard kill:

- `S_0-D_0^C` was positive in approximately 96.5% and negative in approximately 3.5% of draws;
- `S_1-D_1^C` was positive in approximately 83.7% and negative in approximately 16.3% of draws.

Thus both decentralized under-provision and over-provision are possible.

A fixed-cost stress test, drawing `F` from `[0,1.25*max{|D_0|,|D_1|,|S_0|,|S_1|}]`, found examples of all major mismatch directions, including:

- no-entry equilibrium while the national planner chooses one gateway;
- no-entry equilibrium while the planner chooses two gateways;
- full-entry equilibrium while the planner chooses one gateway;
- full-entry equilibrium while the planner chooses no gateways;
- asymmetric one-gateway equilibrium while the planner chooses two;
- asymmetric one-gateway equilibrium while the planner chooses zero.

These are preview results only; no welfare theorem is elevated before novelty re-kill.

## Two simple benchmark calibrations

At `A=1, kappa=5, rho=1.5, tau=1.2, beta=0.5`:

- `D_0^C ≈ 0.056763`,
- `D_1^C ≈ 0.045735`,
- `Gamma^C ≈ -0.011028`,
- `S_0 ≈ 0.077070`,
- `S_1 ≈ 0.055015`.

At the same values but `beta=1.2`:

- `D_0^C ≈ 0.475765`,
- `D_1^C ≈ 0.538184`,
- `Gamma^C ≈ +0.062419`,
- `S_0 ≈ 0.862707`,
- `S_1 ≈ 0.987544`.

The first is gateway substitution; the second is gateway complementarity.