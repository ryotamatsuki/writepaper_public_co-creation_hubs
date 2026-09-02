# Frozen Model Contract

Input main SHA: `3e525ecedfaa4fbac56652d7a40c3f6adf187d7e`.

The Stage 3 CSDS architecture is preserved:

- two regions and two local governments;
- one public matching Hub per region;
- one metropolitan private intermediary;
- two unit masses of immobile projects, one in each region;
- `theta ~ U[0,1]` in each region;
- `w_1(theta)=c+u(1-theta)`, `w_2(theta)=c+u theta`;
- `M_i=x_i w_i(theta)`, `x_i in [0,1]`;
- fixed metro project-side match value `m`;
- common real use cost `kappa`;
- `U_i=M_i-kappa`, `U_T=m-p_T-kappa`, `U_0=0`;
- public investment cost `gamma x_i^2/2`;
- metro profit `p_T D_T`, marginal monetary cost normalized to zero;
- timing `(x_1,x_2) -> p_T -> project choice -> matching`;
- SPNE.

No dynamics, referral, Cournot, taste shocks, relocation, public prices, endogenous metro quality, bargaining, partner entry, or strategic ecosystem design is added.

Stage 3 historical files are not rewritten. Stage 4 records corrections or failures separately.