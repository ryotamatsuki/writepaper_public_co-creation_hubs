# Minimal Model

Parameters:

- `p in (0,1)`: partner-specific compatibility probability.
- `0<l<m<h<a`: low/common, rival-specialist and preferred-specialist productivity gains; `a=A-c` is residual market size.
- `sigma>=0`: resident-partner share of the same match value, normalized as `sigma*x` for a match producing beneficiary productivity `x`.
- `F>0`: binary public hub cost.

Each region has one immobile beneficiary firm and one project. Project type is independently `theta_1` or `theta_2` with probability one half.

The matching stage is the finite R2 architecture documented in `PARTNER_SET_FORMALIZATION.md`.

A match with productivity gain `x` lowers beneficiary marginal cost from `c` to `c-x`. Partner-side collaboration surplus is `sigma*x` and accrues to the partner's region.

Downstream demand is `P=A-Q`. The condition `a>h` guarantees positive Cournot quantities in every matching state.

No price, Hotelling term, network utility, capacity, multihoming, relocation, bargaining or dynamic state is used.
