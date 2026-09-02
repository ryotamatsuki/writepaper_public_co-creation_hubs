# Nested Benchmark Implementation

All benchmarks use the same `theta,c,u,m,kappa,gamma` primitives.

## B1 — one local public + metro private
Remove `G_j/H_j`; retain `x_i`, `p_T`, `w_i`, metro `m`.

## B2 — two local public, no metro
Remove `H_T`; firms choose `H_1,H_2,0`; public-public threshold is inherited unchanged.

## B3 — two local public + fixed metro outside option
Keep `H_T` but set `p_T=bar p`; no private optimization.

## G — full game
Both `x_1,x_2` and follower `p_T` are endogenous.

These are exact restrictions/player removals of the selected candidate model. They are not claimed to reproduce Takahashi, Kim, Liu or Damiano–Li exactly.