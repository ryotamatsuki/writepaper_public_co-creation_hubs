# Stage 33 Lean verification (proof-critical core)

This branch is separate from the production frozen theory on main. The input model is the Stage 33 two-region public/private innovation-intermediation game.

## Reproducibility

Pinned Lean: v4.32.1. Pinned Mathlib commit: 520045ab14e26149ee970e2e617ca04b09bde5d6.

Run from this directory:

```bash
lake update
lake exe cache get
lake build --wfail
```

The GitHub Actions workflow also rejects `sorry`, `admit`, and custom axioms.

## Scope map

- **Sorting**: tail mass and first moment, positive sorting wedge, global effect of a higher interior ability cutoff, positive-ability comparative statics for interior fees.
- **Welfare**: exact regional-vs-restricted-planner public-entry threshold and positive/negative foreign externality.
- **Entry threshold**: polynomial equivalence for a private-only-to-both public entry decision.
- **Entry equilibrium**: exact Boolean logic for four binary deviation gains.

**Not proved by this file**: Cournot payoff derivation from demand, research Nash/effort, endogenous full-domain private-private pricing equilibrium, global 16-profile SPNE with all contingent price games, social optimality across all configurations, empirical applicability, or newness relative to the literature.

The economic model and full strategic proof still require independent mathematical certification. The proof of an encoded implication does not establish its economic antecedents.
