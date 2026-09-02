# Homogeneous Cournot Gateway — Decision Log

Input: PR #14 head `361e7c9ce1b6ed66a0cc9af545f8dc23228b88cc`.

The Stage 4R participation-payoff gateway model is preserved as a **DIAGNOSTIC PRECURSOR ONLY**. It is not treated as the canonical economic model because beneficiary firms did not compete explicitly in a product market and their profits were not product-market equilibrium profits.

## Canonical replacement

The candidate architecture is now:

`regional government gateway choice -> firm ecosystem engagement -> homogeneous Cournot competition`.

Firms produce perfectly substitutable goods. Gateway provision changes only the cost of ecosystem engagement. Engagement changes marginal cost. Marginal costs determine Cournot quantities, market price, market shares and endogenous profits.

## Stage 3C verdict

**PASS — ENDOGENOUS COMPETITION ARCHITECTURE COHERENT.**

No differentiated-product parameter, firm entry, extra spillover, capacity constraint, subsidy or additional strategic stage is required.

## Stage 4C verdict

**GO — HOMOGENEOUS COURNOT PRESERVES ENDOGENOUS GATEWAY STRATEGIC REVERSAL.**

Using the regional-welfare baseline `W_i = pi_i + CS/2 - F E_i`, the exact deletion result is:

- at `beta=0`, `Gamma^C < 0` throughout the canonical stable domain;
- an exact interior calibration with `beta=6/5` yields `Gamma^C>0`;
- an exact interior calibration with `beta=1/2` yields `Gamma^C<0`;
- strict local perturbations preserve both signs, so both are nonempty open sets.

Therefore product differentiation is **NOT NECESSARY**.

## Objective comparison

Producer-only regional government `W_i^P=pi_i-FE_i` also admits both signs. However, producer-only incentives can produce `Gamma^C>0` even at `beta=0` because of strategic rent shifting/business stealing. Regional welfare with half of homogeneous-market consumer surplus removes that clean identification problem: in that baseline, network thickening is necessary for gateway complementarity.

Classification: **RESULT CHANGES QUALITATIVELY** across government objectives.

## Prior model classification

The Stage 4R model is **PARTIAL DIAGNOSTIC ONLY**. It correctly isolated gateway overlap and common-network effects, but product-market rivalry creates a genuine additional force and eliminates the old baseline-core-attractiveness threshold: in the Cournot model `A=a-cbar` scales gross payoffs by `A^2` and does not determine the sign of `Gamma^C`.

## Routing

Freeze the homogeneous-Cournot gateway model as the canonical candidate and send this exact model—not differentiated Cournot—to proposition-level prior-art re-kill.