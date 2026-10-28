# Stage 33 — Lean 4 proof-critical formal verification certificate

Date: 2026-10-28 (Japan)
Repository: `ryotamatsuki/writepaper_public_co-creation_hubs`
Research base: `research/hotelling-entry-experimentation`
Formal source branch: `formal/stage33-lean-verification`
Validated formal source commit: `61d6ecf2fa2fdf70428518d099f24f531e715b81`
Pull request: https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/pull/50

## Formal verification state: PASS — SELECTED PROOF-CRITICAL CORE ONLY

The pinned Lean build and no-hole audit completed successfully:

- GitHub Actions workflow: `stage33-lean-formal`.
- Passing push run: https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/40609975310
- Run ID: `40609975310`; job ID: `122558868954`.
- Build result: `Build completed successfully (8656 jobs).`
- Lean: v4.32.1.
- Mathlib commit: `520045ab14e26149ee970e2e617ca04b09bde5d6`.
- No project `sorry`, `admit`, or custom `axiom`: `LEAN_ESCAPE_HATCH_AUDIT=PASS`.
- `#print axioms` on core claims prints only `propext`, `Classical.choice`, `Quot.sound`, i.e. standard Lean/Mathlib logical axioms.
- The CI was actually executed; this is **not** an uncompiled Lean draft, Z3-only evidence, or an anticipated result.

## Claim to Lean theorem mapping

| Paper component | Lean theorem(s) | Formalized implication | Not certified by this component |
|---|---|---|---|
| Uniform participant skill composition | `skill_factorization`, `selection_wedge`, `selection_wedge_strict`, `new_participant_skill` | Exact tail moments and the positive skill premium of cutoff sorting over random selection | Derivation of that cutoff from all participation subgames |
| Conditional sorting comparisons | `private_fee_increases`, `cutoff_increases`, `participation_decreases`, `aggregate_skill_decreases` | Positive-U, positive-T interior-fee/cutoff implications and decrease of mass/first moment as cutoff rises | Global optimal pricing in all demand regimes, arbitrary ability distributions and two-private-hub fee equilibrium |
| Local/planner welfare wedge | `planner_local_wedge`, `positive_foreign_externality`, `negative_foreign_externality` | Exact public-entry foreign spillover and both possible welfare-wedge signs | Correctness of underlying downstream Cournot payoff primitives or all 16-profile planner global optimum |
| Conditional public-entry threshold | `public_entry_polynomial_identity`, `public_entry_iff_polynomial` | Exact equivalence between `F_G <= local public-entry value` and nonnegative cubic entry polynomial for T>0 | Proving the actual equilibrium continuation for every entry deviation |
| Baseline interval signs | `cubic_positive_on_half`, `baseline_phi_positive`, `baseline_entry_negative_at_zero`, `baseline_entry_positive_at_half` | For the exact rational Stage 33 coefficient inputs, Phi is positive for every lambda in [0,1/2] and Psi crosses signs at interval endpoints | Derivation of these coefficients directly from Cournot primitives *inside Lean*, formal differentiation/IVT/root uniqueness from the model |
| Binary entry logic | `four_player_candidate_NE_iff` | With four given one-flip gains, all nonpositive iff no profitable one-player binary flip | Whether those gains are the correct downstream payoffs, exclusion of all 14 other profiles, mixed equilibrium |

**The theorem labels and numeric coefficient inputs are mapped to the existing Stage 33 derivations and exact-rational Python certificate; that mapping was not formalized independently from market primitives.** The formal verification proves these encoded implications, **not** the entire five-stage game.

## Baseline exact-rational scope

The coefficient inputs for `phi0`–`phi3` and `psi0`–`psi3` are exact rationals copied from `stage8_formal/meta/stage33_exact_rational_certificate.json` of the preceding SMT/exact-rational package (model stage33, dated 2026-10-28). Lean `norm_num` checks their signs and required endpoint inequalities; no floating-point approximations are used inside the Lean proof. The source-to-coefficient provenance remains an explicitly separate audit.

Strict positivity of Phi is formalized on the closed interval [0,1/2]; full uniqueness of the Stage 33 entry threshold follows with additional analytic differentiation and the real-analysis intermediate value argument, not with a complete Lean reconstruction of the entire economic objective.

## Missing from end-to-end formal coverage

1. Closed-form Cournot equilibrium and JV effort derived from primitive demands and the full admissible strategy domains.
2. Existence, global best responses, possible ties and multiplicity in the two-private-hub endogenous fee game outside the Stage 4A sufficient-condition region.
3. All equilibrium continuations and 16-profile global SPNE/first-stage equilibrium correspondence.
4. Comparison against general alternative ability distributions and other product-market regimes.
5. The economy-to-formal-definition mapping of all primitives, and the economic novelty/institutional claims.
6. Independent author certification of proof and economic interpretation.

Thus the strongest valid formal claim is `FORMAL VERIFICATION PASS — PROOF-CRITICAL CORE`, **not** `FULL GAME FORMALLY VERIFIED`, `GLOBAL PRICE EQUILIBRIUM VERIFIED IN LEAN`, or `THEORY FROZEN`.

## Reproduction

At the above source commit:

```bash
cd formal
lake update
lake exe cache get
lake build --wfail
```

CI source: `.github/workflows/stage33-lean.yml`. Formal source: `formal/Stage33Formal.lean`. The CI also greps for custom axioms and proof holes. The pinned build log and axiom prints can be checked in the run linked above.

## Stage 8 / author responsibility

This certificate fulfills the targeted Lean verification operation. It does not constitute the author's intellectual-contribution/sign-off record. Stage 8 final theory freeze remains subject to author review of claims, proof scope and contributions, preservation of independent Stage 4A/Stage 7.5A evidence, and deliberate acceptance of the frozen manuscript's scope. No branch has been merged into main.
