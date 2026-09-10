# Partner Surplus Derivation

The Stage-3/4 architecture has partner cost `c~U[0,1]` and a gross per-platform participation benefit `r_h`. A partner can multi-home, so surplus is additive across platforms.

The participation mass is

`m_h = clip(r_h,0,1)`.

The gross benefit itself is **not** clipped. National partner surplus at platform `h` is therefore

`S_h = ∫_0^{m_h}(r_h-c)dc = r_h m_h - m_h^2/2`.

Equivalently,

- if `r_h <= 0`, then `m_h=0` and `S_h=0`;
- if `0 < r_h < 1`, then `m_h=r_h` and `S_h=r_h^2/2`;
- if `r_h >= 1`, then `m_h=1` and `S_h=r_h-1/2`.

With two symmetric residence groups, each region owns half of the national partner population and therefore receives `S_h/2` from platform `h`.

Hence the familiar regional shortcut `m_h^2/4` is valid only when support participation is interior (`0<r_h<1`, so `m_h=r_h`). It must not be extended to saturation. The all-regime evaluator uses the primitive integral `0.5*(r_h*m_h-m_h^2/2)` for each region and each platform.

The piecewise formula is continuous at both clipping boundaries. At `r_h=0`, both adjacent expressions give zero; at `r_h=1`, both the interior expression and the saturated expression give national surplus `1/2` (regional surplus `1/4`).

This correction does not change partner utility, the support participation cap, or the interior formulas used by the local analytic theorem. It corrects only welfare accounting when an off-path history enters a support-saturation region.
