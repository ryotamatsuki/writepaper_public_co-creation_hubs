# B3 Best-Response Derivation

Fix `p=bar p`. Let `W_i(x_i,x_j,bar p)` denote the participation-reduced regional objective and

`F_i^{B3}=W_{i,x_i}`.

On an interior differentiable best-response branch,

`BR_i^{B3′}= -F_{i,x_j}^{B3}/F_{i,x_i}^{B3} = -W_{i,x_ix_j}/W_{i,x_ix_i}`.

If the explicit B3 SOC condition `W_{i,x_ix_i}<0` holds, then

`sign(BR_i^{B3′})=sign(M_B3)`, where `M_B3:=W_{i,x_ix_j}`.

## Economic block

At fixed price the only cross-region connection in the central regime is the shared two-sided participation system. The chain is

`x_j ↑ -> n_j^P ↑ -> n_j^F ↑ -> n_T^F ↓ -> n_T^P ↓ -> q_T ↓ -> region-i projects shift toward H_i -> marginal value of x_i ↑`.

The exact participation derivatives in `PARTICIPATION_REDUCTION.md` make this channel analytic. In particular, at `beta=0` (`delta=0`) the cross participation terms vanish, recovering the production result `M_B3=0` exactly.

## Limit of the reduction

For `beta>0`, fully eliminating the equilibrium cutoffs from `M_B3` produces a high-degree rational expression. The sign is therefore retained as an exact reduced derivative block unless a transparent primitive factorization can be found. At the canonical B3 equilibrium it is strictly positive, but the earlier counterexample search does not justify declaring it positive throughout the entire primitive domain.

**Gate 2 for B3: PASS in exact reduced-form block representation. Primitive all-regime sign: not proved.**