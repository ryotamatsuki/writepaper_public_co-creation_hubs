# Ecosystem Engagement Subgame

Let

`m=2+beta`, `n=beta-1`, `b=2mn`.

For a gateway profile producing engagement-cost coefficients `(K_i,K_j)`, define

`a_i = 9K_i-2m^2`,

`a_j = 9K_j-2m^2`,

`Delta = a_i a_j-b^2`.

After substituting the Cournot equilibrium into profit, firm `i` solves

`max_{x_i>=0} q_i(x_i,x_j)^2-(K_i/2)x_i^2`.

The interior FOC is

`a_i x_i - b x_j = 2mA`.

The second FOC is symmetric. Hence

`x_i^* = 2mA(a_j+b)/Delta`,

`x_j^* = 2mA(a_i+b)/Delta`.

The corresponding Cournot quantities can be written compactly as

`q_i^* = 3A K_i(a_j+b)/Delta`,

`q_j^* = 3A K_j(a_i+b)/Delta`.

The FOC also implies

`q_i^* = 3K_i x_i^*/(2m)`.

Therefore equilibrium total firm profit is

`pi_i^* = A^2 K_i a_i (a_j+b)^2 / Delta^2`.

## Engagement strategic interaction

Firm `i`'s engagement best response is

`x_i = [2mA+b x_j]/a_i`.

Thus its slope is

`b/a_i`.

Under the stable domain `a_i>0`, the sign is the sign of `beta-1`:

- `beta<1`: firm engagement choices are strategic substitutes;
- `beta>1`: firm engagement choices are strategic complements;
- `beta=1`: engagement best responses are locally independent.

This is distinct from the local-government gateway interaction.

## Gateway profiles

Use

`H=kappa`,

`L=kappa-rho`,

`M=kappa-rho+tau`.

Then the four firm subgames are obtained from the same closed form by evaluating:

- `00: (K_i,K_j)=(H,H)`;
- `10: (L,M)`;
- `01: (M,L)`;
- `11: (L,L)`.

No separate solution concept or additional state variable is needed.