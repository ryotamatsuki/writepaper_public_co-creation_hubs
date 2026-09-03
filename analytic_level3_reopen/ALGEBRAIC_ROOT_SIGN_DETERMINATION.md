# Algebraic-Root Sign Determination

L3-1 closes the two canonical sign obligations that remained open in the first PR #40 draft.

The proof does not solve the degree-31 root by radicals. `code/verify_l31.py` first isolates the coupled G/B3 stationary solution with an exact rational Krawczyk certificate. It then computes reduced Hessians through the participation/private-price implicit systems using exact `Fraction` interval arithmetic and cofactor/determinant inverses.

Certified intervals:

- `H11_B3 in (-0.118357218885676,-0.118357218885651)`;
- `H12_B3 in (0.0135431333332136,0.0135431333332153)`;
- `H11_G in (-0.199391062101166,-0.199391062101008)`;
- `H12_G in (-0.00480853615115529,-0.00480853615112228)`;
- private price SOC in `(-52.1627947003595,-52.1627947003591)`;
- `p_{x_i}` in `(-0.0165046011641228,-0.0165046011641190)`.

Hence, under strict SOCs, the BR-slope intervals are approximately B3 `(0.114425917242068,0.114425917242105)` and G `(-0.024116106812849,-0.024116106812664)`.

This is an exact sign-at-algebraic-root result for the canonical calibration. What remains missing for Level 3 is the non-calibrated primitive projection/root-selection theorem.
