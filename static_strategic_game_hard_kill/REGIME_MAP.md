# Regime Map

The exact regime map can be represented by the public upper envelope `B(theta)=max{V_1(theta),V_2(theta)}` and a constant middle alternative.

Let `r=q` when `q>kappa`; otherwise let `r=kappa` and metro is inactive.

With an interior public crossing:

1. **Public-public split:** `r<=K`. The middle constant option never wins.
2. **Public–middle–public:** `K<r<min{x_1d,x_2d}`. Both public tails survive and the constant option occupies a middle interval.
3. **One public tail + middle:** `r` exceeds one public endpoint but not the other.
4. **Middle monopoly:** `r>=max{x_1d,x_2d}`.

When `r=q>kappa`, “middle” means metro. When `r=kappa>=q`, “middle” means non-participation.

If the public crossing is outside `[0,1]`, one public Hub dominates the other everywhere and the game reduces locally to a one-public-vs-constant problem.

This clipped-envelope representation covers the requested all-public, public-public, public-metro-public, public-metro, metro-only, public-plus-zero, and zero-use corner cases without pretending every verbal label is a separate economic regime.

Private profit is continuous and piecewise quadratic in `q`; the global private optimum is obtained by checking the stationary point on each clipping segment and the kink points. Public objectives are continuous but piecewise differentiable. The central interior regime is therefore analytically solvable; the global game is algorithmically tractable rather than represented by one global polynomial.