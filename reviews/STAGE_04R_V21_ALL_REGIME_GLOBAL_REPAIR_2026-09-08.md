# Stage 4R — v2.1 All-Regime Global Equilibrium Repair

Date: 2026-09-08 JST

Repository: `ryotamatsuki/writepaper_public_co-creation_hubs`

Workflow authority: `ryotamatsuki/research-paper-workflow` v2.1 @ `fe6fa5f2a94632be578d5e9fab39e7df86013113`

Template: `templates/STAGE_04_MINIMAL_MODEL.md`

Predecessor: Stage 4A v2.0 NO-GO / REOPEN STAGE 4, which found a profitable finite off-regime public deviation at the old numerical witness `(beta,gamma,tau)=(0.05,0.9,0.05)`.

## 1. Executive verdict

**GO** — construction-level only. Route immediately to Stage 4A independent mathematical adversarial certification.

The model itself is unchanged. The repair re-solves the stated strategy game across all primitive project routes and public finite deviations, then searches the existing primitive parameter space for a region where the local complements-to-substitutes mechanism is also supported by global public best responses.

A construction-level passing witness is:

- `beta = 0.01`;
- `gamma = 0.825`;
- `tau = 0.35`;
- all other primitives unchanged from the existing model.

The larger `tau` is not used as an ad hoc new friction: `tau` was already a primitive. At `tau=.35`, `kappa_L+tau=.62 > v+alpha=.60`, so the nonresident rival public route is globally dominated by the outside option at every history. This removes the specific cross-region free-riding regime that killed the old witness. The lower network-effect strength and lower cost curvature are likewise values of existing primitives, not model changes.

## 2. Repaired symmetric candidates

### Full game G

Construction-stage stationary candidate:

- `x_G ~= 0.8371022382`;
- locally computed private price `p_G ~= 0.0184036834`;
- all-regime globally re-solved private price `p_G ~= 0.0184036832`;
- public own second derivative `H_11 ~= -0.21031634 < 0`;
- public cross derivative `H_12 ~= -0.00417719 < 0`;
- best-response slope `BR_G' ~= -0.01986145 < 0`.

Dense all-regime public best-response search, holding the rival at the candidate and re-solving downstream private pricing at every public deviation, returns

- `BR_G(x_G) ~= 0.8370893604`;
- candidate welfare `~= 0.29487864055`;
- best detected welfare `~= 0.29487864056`;
- gain over candidate `~= 9.4e-12`.

Within construction tolerance, the candidate is the global public best response.

### Matched-price B3

With the private fee fixed at the all-regime G price:

- `x_B3 ~= 0.8258903860`;
- public own second derivative `H_11 ~= -0.18918949 < 0`;
- public cross derivative `H_12 ~= +0.000750729 > 0`;
- best-response slope `BR_B3' ~= +0.00396813 > 0`.

Dense all-regime public best-response search returns

- `BR_B3(x_B3) ~= 0.8258866749`;
- candidate welfare `~= 0.28933428339`;
- best detected welfare `~= 0.28933428339`;
- reported gain over candidate is numerical noise (`~-1.3e-12`).

Thus the repaired construction witness preserves

`BR_G' < 0 < BR_B3'`

while also surviving finite/global public deviation search.

## 3. Negative repair diagnostics retained

The repair search did not cherry-pick the first convenient parameter change.

- `tau` alone failed: even after making the rival public route globally dominated, low-investment deviations to the private route remained profitable at the old `(beta,gamma)`.
- `gamma` alone failed over the searched range.
- `rho` alone failed over the searched range.
- At `gamma=.85`, small beta values preserved the local sign reversal but B3 still had a profitable low-investment deviation.
- The `(beta,gamma)=(.01,.825)` witness is selected from a broader screen in which several `gamma=.825` small-beta points passed the construction-level global search.

These failed repairs remain regression evidence and are not deleted.

## 4. Continuation semantics

The all-regime evaluator:

- reconstructs project choice from the primitive upper envelope over `0,H1,H2,HT`;
- applies partner-participation clipping from `c~U[0,1]`;
- solves participation directly;
- searches the economically relevant private-price interval globally;
- never maps `None`, branch failure, nonconvergence, or invalid active sets to low profit;
- treats unresolved/multiple material continuations as blockers.

The selected candidate did not encounter a material unresolved continuation in the dense construction check.

## 5. Claim status after repair

Construction-level status only:

| Claim | Stage 4R status |
|---|---|
| Existing model admits a regular local B3-complements/G-substitutes pair | PASS |
| Same pair survives all-regime finite/global public deviation search at selected witness | PASS, numerical construction evidence |
| Private continuation at selected on-path G candidate is globally solved numerically | PASS, construction evidence |
| Global Nash/SPNE status as a theorem/certified object | NOT YET CERTIFIED — Stage 4A required |
| Old exact L3-1 certificate at `(beta,gamma,tau)=(.05,.9,.05)` | STALE for the repaired witness |
| Old Stage 8–14 frozen package | STALE pending new Stage 4A and downstream recertification |

## 6. Canonical Stage 4 verdict

`GO`

## 7. Route

Proceed to **Stage 4A — Independent Mathematical Adversarial Certification Gate**.

Stage 4A must not import or call the Stage 4 repair solver. It must independently reconstruct the selected witness, attack the complete public strategy domain, verify downstream continuation semantics, and attempt to falsify both global equilibrium status and the sign-reversal claim.
