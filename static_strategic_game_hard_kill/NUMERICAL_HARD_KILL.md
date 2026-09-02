# Numerical Hard-Kill

Reproducible seed: `20260902`. Python 3.13.5. Sign tolerance: `1e-10`.

The saved `code/numerical_hard_kill.py` reproduces both experiments below exactly.

## Experiment A — 100,000 admissible central-regime states

Rejection sampling continued until 100,000 states satisfied the interior public crossing and public–metro–public validity conditions. Under the saved sampling domain, 158,039 raw draws were required.

Counts:

- `dp_T*/dx_1 < 0`: 100,000
- G public cross-partial `<0`: 99,767
- G public cross-partial `>0`: 233
- B2 public cross-partial `<0`: 98,835
- B2 public cross-partial `>0`: 1,165
- B2 complement / G substitute: 932
- B2 substitute / G complement: 0

Thus off-symmetric complementarity is possible in G, but this broad state search did not reveal G creating complementarity where B2 was substitutive.

## Experiment B — 100,000 constructed valid symmetric G stationary equilibria

For each draw, `c,u,x,m,kappa` were sampled strictly inside the active symmetric public–metro–public domain and `gamma` was backed out from the exact G symmetric Nash FOC. The G SOC was then checked. All 100,000 raw draws were accepted.

At every draw:

- the fixed-offer B3 marginal evaluated at the G stationary equilibrium was positive: 100,000/100,000;
- hence removing endogenous private price response creates upward local pressure on public investment.

Under the strict frozen social-welfare accounting, the planner marginal at the G Nash point was:

- positive in 7,250 cases;
- negative in 92,750 cases.

The sign classification matches the analytical `y*` threshold in `DECENTRALIZATION_WEDGE.md`.

For the subset where both B2 Nash and planner formulas were interior, inside `[0,1]`, and full participation held (48,140 draws):

- B2 overinvestment: 39,792;
- B2 underinvestment: 8,348.

The experiments are diagnostic only. Exact propositions and the decisive characteristic-space isomorphism are symbolic results, not numerical inferences.