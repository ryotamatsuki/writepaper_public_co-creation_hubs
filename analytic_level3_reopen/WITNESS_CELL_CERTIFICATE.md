# Canonical Witness Certificate

The canonical finite-decimal primitives are represented exactly as rationals. L3-1 now certifies the **coupled** G/B3 stationary root, not G alone.

An exact rational Krawczyk operator proves one and only one solution inside:

- `t_G in (0.5848900730,0.5848900731)`;
- `q_G in (0.1388515735,0.1388515736)`;
- `p_G in (0.0227466514,0.0227466515)`;
- `t_B in (0.6100908625,0.6100908626)`;
- `q_B in (0.1402669351,0.1402669352)`.

Reconstruction gives `x_G≈0.684028196361275` and `x_B≈0.656020390747393`, so the matched-price coupling is respected with `x_B != x_G` and common algebraic `p_G`.

Exact rational interval implicit differentiation then proves:

- B3 public SOC `<0` and cross derivative `>0`;
- G public SOC `<0` and cross derivative `<0`;
- private price SOC `<0`;
- `p_{x_i}<0`;
- follower/public Jacobians nonsingular on the isolating box.

The canonical reversal is therefore an exact computer-assisted algebraic certificate, not a floating-point witness.

This is **not** a primitive-parameter witness-cell theorem: a full CAD/sign-invariant primitive cell around the calibration was not computed. That distinction prevents overclaiming Level 3.
