# Stage 4R-MC-G — Executive Verdict

## Verdict

**NO-GO — MATCHING STRUCTURE COLLAPSES**

The CSDS game is analytically tractable and the full game does generate a genuine price-mediated public-public cross-effect. However, the selected matching primitive fails the Stage 4 matching-specificity gate.

For `d=c+u`, define `Q_i=x_i d` and `t_i=x_i u`. Then, for `theta in [0,1]`,

- `U_1(theta)=Q_1-t_1 theta-kappa`,
- `U_2(theta)=Q_2-t_2(1-theta)-kappa`,
- `U_T=m-p_T-kappa`.

Thus the complete project-choice subgame is exactly a one-dimensional linear characteristic-space / Hotelling-type facility-choice game with provider-specific slopes constrained by `t_i=(u/d)Q_i`. Risk-neutral expected matching enters only through deterministic expected utility. No stochastic matching state, search sequence, partner draw, or other matching-specific strategic object remains.

## Diagnostic results that survive algebraically

In the interior public–metro–public regime, let `a=2c+u` and

`K=a x_1 x_2/(x_1+x_2)`.

The private best response is

`p_T*=(m-K)/2`, `q_T*=m-p_T*=(m+K)/2`,

with `dp_T*/dx_i<0`.

With a fixed metro offer (B3), the regional-government cross-partial is exactly zero. With endogenous metro pricing (G),

`W_{1,12}^G = -a[a x_1 x_2(x_1+7x_2)+2m(x_2^2-x_1^2)]/[4u(x_1+x_2)^4]`.

At symmetry this is `-a^2/(8ux)<0`. Hence endogenous metro pricing re-couples the two governments. This is a real full-game feedback, but it is generated entirely by the characteristic-space payoff geometry and therefore does not rescue the matching contribution.

## Additionality failure in the active metro regime

If the metro intermediary is active, `q_T*>kappa`; because its payoff is type-independent, every project not choosing a public Hub strictly prefers metro to non-participation. Hence the active public–metro–public regime has full participation and marginal local-Hub expansion is metro displacement, not `0 -> H_i` market creation. The model cannot simultaneously deliver active metro competition and the extensive additionality margin that motivated the research question.

## Routing

Do not proceed to Stage 5/6 for this static CSDS model. Do not activate the reserved dynamic ecosystem extension automatically. A separate diagnosis is required before any dynamic reboot.