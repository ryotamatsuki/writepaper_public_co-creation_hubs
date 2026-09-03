# Parameter Restriction Register

Baseline: `v>0`, `alpha>0`, `beta>=0`, `rho>=0`, `rho_T>=0`, `gamma>0`, `kappa_L>0`, `kappa_T>=0`, `tau>=0`; `x_i in [0,1]`, `p_T>=0`.

Interior partner block requires `rho+x_i+beta n_i^F in (0,1)` and `rho_T+beta n_T^F in (0,1)`.

Central smooth allocation regime for each home region requires:
- `0<s<t_i<1`;
- `q_i>q_T`;
- `kappa_T+p_T<kappa_L`;
- `s=(kappa_T+p_T)/q_T`;
- `t_i=(kappa_L-kappa_T-p_T)/(q_i-q_T)`.

The certified equilibrium requires negative private and public own SOCs and nonsingular fixed-point/equilibrium Jacobians. The local IFT robustness claims additionally require positive density at relevant cutoffs and persistence of regularity/SOCs.

P2-R4 is defined by strict sign inequalities `M_B3>0` and `M_B3+P_price<0`. No simpler global parameter restriction is frozen.