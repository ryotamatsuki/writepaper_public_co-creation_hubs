# Stage 4A v2.0 — Retrospective Independent Mathematical Certification

Date: 2026-09-07

Workflow authority: `ryotamatsuki/research-paper-workflow` current `main`, v2.0 architecture, commit `698b6a42cc177b30300005c1b4e25956b2cf79d9`.

Audit base: Stage-14 QA branch head `58f13f23bc18989a6404325e21068b339c45a057`. Stage-14 changes are submission-only; the active mathematical content is the merged Stage-13 theorem/certificate hierarchy.

## Independence contract

`code/independent_adversarial_audit.py` does not import or call the production central-regime solver, `numerical_hard_kill.py`, or the L3-1 Krawczyk verifier. It reconstructs the economic allocation directly from primitive utilities:

1. route demand is obtained from the upper envelope of `0,H1,H2,HT` affine project utilities on `z in [0,1]`;
2. partner participation is reconstructed from the primitive `c~U[0,1]` rule, with clipping to `[0,1]`;
3. the participation fixed point is solved directly from the induced demand mapping and checked from multiple starts at the audited states;
4. the private price is searched on the full economically relevant interval. Since `n_T^P<=1`, `q_T<=v+alpha=0.6`; therefore `p_T>=0.58` implies zero private demand at the canonical primitives, so any positive-profit global optimum lies in `[0,0.58]`;
5. regional welfare is integrated directly from the realized upper envelope, plus primitive partner surplus and public cost.

The same file separately re-derives the beta-zero analytic identities with SymPy from the primitive central-regime equations rather than copying the production expressions.

## Adversarial result

The local analytic identities survive independent reconstruction, but the canonical G and B3 stationary configurations do not survive finite-deviation/global-best-response attack.

At the canonical full-game stationary configuration:

- `x_1=x_2 ~= 0.684028196361275`;
- independently re-solved private price `p_T ~= 0.02274665143`;
- direct regional welfare `W_1 ~= 0.22301946658`.

Holding `x_2` fixed, government 1 can deviate to `x_1=0.1875`. Re-solving the private price and participation globally gives approximately:

- `p_T ~= 0.02518788380`;
- project masses `(n_1^F,n_2^F,n_T^F) ~= (0,0.734535,0.622508)`;
- `W_1 ~= 0.25444926939`;
- gain `~= +0.03142980281`.

The deviation is economically transparent. At `x_1=0.1875` and `n_1^F=0`,

`q_1 = v + alpha(rho+x_1) = 0.26875 < kappa_L = 0.27`,

so the deviating region's own public hub is below the outside option even for `z=1`; high-type projects instead use the rival public hub. In the own-hub-inactive interval, the only direct `x_1` term is own partner surplus minus public cost, whose interior maximizer is exactly

`x_1 = rho/(2 gamma - 1) = 0.15/0.8 = 0.1875`.

The matched-price B3 stationary configuration fails the same test:

- candidate `W_1 ~= 0.21091151931` at `x_1=x_2 ~= 0.656020390747393`;
- deviation `x_1=0.1875` gives `W_1 ~= 0.23886596269`;
- gain `~= +0.02795444338`.

## Consequence

The Stage-13 local sign-reversal theorem can remain mathematically valid as a conditional local/regular-branch statement. What fails is the promotion of the canonical stationary pair to public Nash / full-game SPNE equilibrium and any welfare statement whose premise is that canonical pair being a decentralized equilibrium.

Under v2.0 fail-closed semantics, this is a Stage-4 equilibrium/globality failure. Downstream Stage 8--14 states are stale until Stage 4 is repaired/re-solved and Stage 4A is rerun.
