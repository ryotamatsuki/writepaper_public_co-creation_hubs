# Project Choice Cutoffs

Write `q_h=v+alpha n_h^P` and `a_rh=kappa_rh+1{h=T}p_T`.

Outside cutoff: `z_h0^r=a_rh/q_h`.

Pairwise cutoff: `z_hg^r=(a_rh-a_rg)/(q_h-q_g)` whenever `q_h != q_g`.

Main regime for a resident of region i:
- `0` for `z<s`,
- `H_T` for `s<z<t_i`,
- local `H_i` for `t_i<z<1`,
with `s=(kappa_T+p_T)/q_T` and `t_i=(kappa_L-kappa_T-p_T)/(q_i-q_T)`.

These formulas were independently checked in the verification scripts.