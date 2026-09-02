# Parameter Domain

Baseline restrictions are

`0<c<c+u<=1`, `u>0`, `0<m<=1`, `kappa>0`, `gamma>0`, `x_i in [0,1]`, `p_T>=0`.

For compactness define

- `d := c+u`,
- `a := 2c+u`,
- `K(x_1,x_2) := a x_1 x_2/(x_1+x_2)`.

When the public-public crossing lies in `(0,1)`, `K` is the minimum gross payoff on the upper envelope of the two public-Hub match values.

The crossing is interior iff

`x_1 d > x_2 c` and `x_2 d > x_1 c`.

For the central public–metro–public regime with zero private marginal monetary cost, the private best response is valid when

1. `m>K` (positive private price and metro entry),
2. `q*=(m+K)/2>kappa` (metro dominates non-participation),
3. `q*<min{x_1 d,x_2 d}` (both public tails remain active).

Because `K` lies between each local line's endpoint values when the public crossing is interior, `q*>K` automatically places both public-metro cutoffs inside the relevant public branches.

At symmetry `x_1=x_2=x`, these conditions reduce to

`a x/2 < m < x(a/2+u)` and `kappa < (m+a x/2)/2`.

All claimed interior results are restricted to these inequalities. KKT/boundary outcomes are not silently treated as interior.