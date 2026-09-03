# Full-Game Best-Response Derivation

Let `W_i(x_i,x_j,p)` be the participation-reduced regional objective with price treated as an argument, and let `p=p^*(x_i,x_j)` solve the private FOC. Define

`W~_i(x_i,x_j)=W_i(x_i,x_j,p^*(x_i,x_j))`.

The full public FOC is `F_i^G=W~_{i,x_i}=0`. On an interior differentiable best-response branch,

`BR_i^{G′}= -W~_{i,x_ix_j}/W~_{i,x_ix_i}`.

## Exact cross derivative

A complete chain rule gives

`W~_{i,x_ix_j}`

`= W_{i,x_ix_j}`

`+ W_{i,x_ip} p_j`

`+ W_{i,x_jp} p_i`

`+ W_{i,pp} p_i p_j`

`+ W_{i,p} p_{ij}`,

where `p_i=p^*_{x_i}`, `p_j=p^*_{x_j}`, and `p_{ij}=p^*_{x_ix_j}`.

Thus the exact decomposition is

`M_G = M_fixed + P_price`,

with

`M_fixed := W_{i,x_ix_j}`

and

`P_price := W_{i,x_ip}p_j + W_{i,x_jp}p_i + W_{i,pp}p_ip_j + W_{i,p}p_{ij}`.

The production Stage-4 description that the full cross effect contains direct/network terms plus first and second derivatives of the price policy is therefore exact. A one-term price-response story would be incomplete.

If `W~_{i,x_ix_i}<0`, then `sign(BR_i^{G′})=sign(M_G)`.

## Mechanism versus headline comparison

At a fixed `(x_i,x_j,p)` the condition `P_price<-M_fixed` is exactly the condition that endogenous repricing turns a positive fixed-price cross derivative negative **at that same state**.

The headline B3/G result is stronger and different: B3 fixes `bar p=p^G` and then solves for a distinct public equilibrium `x^{B3}`, whereas G is evaluated at `x^G`. Hence `P_price<-M_fixed` at `x^G` is an exact mechanism diagnostic but not, by itself, a necessary-and-sufficient primitive condition for `BR^{B3′}(x^{B3})>0>BR^{G′}(x^G)`.

This distinction is a central L3-0 finding.