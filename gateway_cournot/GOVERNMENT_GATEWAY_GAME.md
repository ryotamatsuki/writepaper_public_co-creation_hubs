# Government Gateway Game

For any pair of engagement-cost coefficients `(K_i,K_j)`, let

`Pi_i(K_i,K_j)=A^2 K_i a_i(a_j+b)^2/Delta^2`,

where

`a_i=9K_i-2(2+beta)^2`,

`a_j=9K_j-2(2+beta)^2`,

`b=2(2+beta)(beta-1)`,

`Delta=a_i a_j-b^2`.

Total output is

`Q(K_i,K_j)=3A[K_i(a_j+b)+K_j(a_i+b)]/Delta`,

so total consumer surplus is

`CS(K_i,K_j)=Q(K_i,K_j)^2/2`.

The canonical regional gross welfare is

`R_i(K_i,K_j)=Pi_i(K_i,K_j)+CS(K_i,K_j)/2`.

Use `H=kappa`, `L=kappa-rho`, `M=kappa-rho+tau`.

Then

- `R_i^{00}=R_i(H,H)`;
- `R_i^{10}=R_i(L,M)`;
- `R_i^{01}=R_i(M,L)`;
- `R_i^{11}=R_i(L,L)`.

Define gross gateway incentives

`D_0^C = R_i(L,M)-R_i(H,H)`,

`D_1^C = R_i(L,L)-R_i(M,L)`,

and

`Gamma^C = D_1^C-D_0^C`.

The exact closed form is a rational function with a high-order numerator. No economically meaningful global three-term factorization was found. The canonical result therefore uses exact compact value functions plus mechanism-deletion proofs rather than presenting a large polynomial as an economic decomposition.

## Gateway-game topology

For fixed gateway cost `F`:

- if `D_1^C<D_0^C`, gateway choices are strategic substitutes. `F<D_1^C` gives full entry, `D_1^C<F<D_0^C` gives asymmetric one-gateway equilibria, and `F>D_0^C` gives no entry (subject to positive-threshold qualifications);
- if `D_0^C<D_1^C`, gateway choices are strategic complements. `F<D_0^C` gives full entry, `D_0^C<F<D_1^C` gives `(0,0)` / `(1,1)` coordination multiplicity, and `F>D_1^C` gives no entry.

The topology is an implication, not a novelty claim.

## Baseline-market profitability

Every equilibrium engagement and output is proportional to `A`, and every gross government payoff difference is proportional to `A^2`. Therefore

`sign(Gamma^C)`

is independent of `A>0`.

`A` scales gateway willingness to pay relative to fixed cost `F`, but it does not determine strategic complementarity versus substitutability.