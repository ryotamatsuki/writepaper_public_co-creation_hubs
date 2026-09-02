# Stage 4M Verdict

## Executive verdict

Stage 3M:

**PASS — COHERENT MINIMAL PRODUCTIVE-MATCHING MICROFOUNDATION.**

Stage 4M:

**GO — PRODUCTIVE-MATCHING PUBLIC-HUB DUALITY SURVIVES THE MINIMAL MODEL.**

## Central results

### T1 — Matching-only substitution

**CONFIRMED.** With referral disabled (`r=0`):

`Gamma_0 = -Delta p[8A+Delta(11-14p^2+7p^3)]/18 < 0`

throughout `A>Delta>0`, `0<p<1`.

### T2 — Referral-induced complementarity

**CONFIRMED.** `Gamma(k)` is strictly increasing in `k=(1-p)^2 r`. If `Gamma_max>0`, there is a unique interior referral threshold. Exact positive and negative interior witnesses exist and are robust to local perturbation.

### T3 — Competition wedge

**CONFIRMED, but refined.** Downstream competition is not necessary for sign reversal. The nonstrategic benchmark has threshold `k=p`. Cournot regional welfare satisfies

`Gamma_C(k=p)=-7 Delta^2 p^2(p^2-4p+5)/18<0`,

so downstream competition strictly raises the referral strength needed for complementarity.

Classification: **ESSENTIAL TO GOVERNMENT WEDGE**.

### T4 — Decentralization wedge

**CONFIRMED.** A planner uses `S_e=D_e+X_e`; regional governments use `D_e`. The first-hub externality satisfies `X_0>0` globally, while the second-hub externality can change sign. Public decentralization therefore changes provision incentives.

## Strategic-sign result

**BOTH SIGNS on nonempty interior open sets.**

Exact calibration `A=1`, `Delta=1/2`, `p=1/10`:

- `r=1/10`: `Gamma=-761087/72000000<0`;
- `r=1/5`: `Gamma=33521/2250000>0`.

Joint ±5% perturbations of all primitives produced 100,000/100,000 sign-preserving valid draws around each witness.

## Large numerical falsification

Seeded 500,000-draw admissible search with `A in [0.2,4]`, `Delta/A in [0.02,0.95]`, `p in [0.005,0.95]`, `r in [0,1]`:

- M0: 0 positive / 500,000 negative / 0 near zero;
- M1: 106,678 positive / 393,322 negative / 0 near zero at tolerance `1e-12`;
- nonstrategic benchmark: 127,377 positive / 372,623 negative;
- 20,699 cases were nonstrategic-complementary but Cournot-substitutive;
- 0 cases were nonstrategic-substitutive but Cournot-complementary.

An independent 10,000-draw enumeration of all four binary matching states verified the closed-form expected regional welfare with maximum absolute error `5.33e-15`.

## Productive-matching primitive verdict

**VALID CANONICAL PRIMITIVE.**

## Shared-pool/referral verdict

**ESSENTIAL** for complementarity in the selected minimal model.

## Immediate absorption result

Caillaud–Jullien, Damiano–Li, Calcagnini et al., and Ichihashi absorb major components, but the immediate one-page specialization test did not recover the full public binary provision + finite productive matching + referral/overlap + downstream competition + regional-welfare threshold package.

This is a Stage 4 survival result only. It is not a novelty claim.

## Previous diagnostic models

- Gateway Stage 4R: **USEFUL DIAGNOSTIC** — identified overlap/access logic but lacked endogenous product-market competition and the institutionally selected matching primitive.
- Cournot Stage 4C: **USEFUL DIAGNOSTIC** — demonstrated the role of endogenous firm competition, but its generic network/R&D-like engagement technology is not the canonical hub primitive after Stage 0M–2M.

Neither old theorem is silently imported.

## Routing

**GO TO STAGE 6M — PROPOSITION-LEVEL NOVELTY RE-KILL.**

Freeze the current minimal model and the exact T1–T4 results. Do not add robustness, differentiation, R&D, dynamics, more regions, or additional matching mechanisms before the exact proposition-level literature re-kill.