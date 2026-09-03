# Canonical Theory Specification

## Research object
Two regional governments choose public facilitation investments while a shared private metropolitan innovation hub subsequently chooses an access fee. The question is whether cross-side network feedback makes public investments complements when the private fee is fixed but substitutes when the private hub optimally reprices.

## Strategic decision-makers and routes
Strategic decision-makers are `G_1`, `G_2`, and the private operator `H_T`. The project route set is `0,H_1,H_2,H_T`. Earlier Stage-4 shorthand that listed `H_1,H_2` among “players” is interpreted as platform/route labels, not additional strategic decision-makers.

## Timing and equilibrium
`(x_1,x_2) -> p_T -> participation fixed point/project sorting -> surplus`. The two governments choose first-stage investments; the private operator chooses `p_T` after observing them. Equilibrium is SPNE.

## Project side
Each region has unit mass `z~U[0,1]`. A focal project single-homes across `0,H_1,H_2,H_T`:
`U_rh(z)=z[v+alpha n_h^P]-kappa_rh-1{h=T}p_T`, `U_r0=0`.
Costs: `kappa_rr=kappa_L`, `kappa_rj=kappa_L+tau`, `kappa_1T=kappa_2T=kappa_T`.
Define `q_h=v+alpha n_h^P`.

## Partner side
Partners have platform-specific cost `c~U[0,1]` and may multihome. In the interior block:
`n_i^P=rho+x_i+beta n_i^F`, `n_T^P=rho_T+beta n_T^F`.

## Controls and costs
`x_i in [0,1]`, public cost `gamma x_i^2/2`. The private operator chooses `p_T>=0` and earns `Pi_T=p_T n_T^F` under zero real private operating cost.

## Smooth central regime
For home region `i`: `0<s<t_i<1`, `q_i>q_T`, `kappa_T+p_T<kappa_L`, with `s=(kappa_T+p_T)/q_T` and `t_i=(kappa_L-kappa_T-p_T)/(q_i-q_T)`. The local ordering is `0 -> H_T -> H_i`.

## Welfare
`PS_i=∫_s^{t_i}(zq_T-kappa_T-p_T)dz+∫_{t_i}^1(zq_i-kappa_L)dz`.
Partner surplus on platform `h` is `b_h^2/2` nationally, `b_h=n_h^P`; under symmetric residence ownership region `i` receives `b_h^2/4`.
`W_i=PS_i+(b_1^2+b_2^2+b_T^2)/4-gamma x_i^2/2`.
`W^N=W_1+W_2+Pi_T`; fee payments cancel as aggregate transfers.

## Benchmarks
B3 retains all three hubs and two-sided participation but fixes `p_T=bar p`; for the headline witness `bar p=p_T^G`. G is the full sequential game with optimized `p_T^*(x_1,x_2)`.

## Headline result
A regular interior witness satisfies `BR_i^{B3′}>0>BR_i^{G′}`. The derivative decomposition is `M_B3>0` and `M_B3+P_price<0`. This is a numerically certified existence result with local continuity/open-set support, not a global sign theorem.

## Essentiality
At `beta=0`, the relevant B3 public-public cross effect is exactly zero. Two-sidedness is therefore essential to the strong complements-to-substitutes comparison, though not to private-price transmission generally. Separate regional objectives are essential to the noncooperative public-public best-response object.

## Welfare content
At the frozen G witness, `dW^N/dx_i≈0.42277>0` after fee-transfer cancellation. This is a real coordination wedge. The strategic-sign reversal is not a welfare-sign reversal.

## Approved robustness
Only local/open-set robustness is approved: deterministic ±0.5% perturbations; small `C^1` distribution perturbations; small `C^1` network-function perturbations; and small asymmetric action-state perturbations under regularity.

## Institutional interpretation
`H_i` is a regional public innovation hub/intermediary; `H_T` is a shared metropolitan private innovation hub; `x_i` is facilitation capacity; `p_T` is an effective paid membership/access price. Project single-homing means one focal project selects one primary route, not firm exclusivity. Partner multihoming is an institutionally plausible abstraction.

## Contribution boundary
The contribution is the fixed-price-complements to endogenous-price-substitutes benchmark reversal. Network complementarity, endogenous price response, two-sidedness, mixed public/private competition, and multiple platforms are known components and are not independently novel.

## Change control
All later manuscript work must conform to this specification. Any substantive alteration requires `THEORY_CHANGE_CONTROL.md` and reopening the affected prior gate.