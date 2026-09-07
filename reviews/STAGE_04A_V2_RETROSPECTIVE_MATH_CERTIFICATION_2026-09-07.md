# Stage 4A — v2.0 Retrospective Independent Mathematical Adversarial Certification

Date: 2026-09-07

Repository: `ryotamatsuki/writepaper_public_co-creation_hubs`

Audit branch: `audit/v2-stage4a-retrospective-certification`

Audit base: `58f13f23bc18989a6404325e21068b339c45a057`

Workflow authority: `ryotamatsuki/research-paper-workflow` current `main`, v2.0 architecture at `698b6a42cc177b30300005c1b4e25956b2cf79d9`.

Template: `templates/STAGE_04A_MATH_RED_TEAM.md`

## 1. Executive adversarial verdict

**NO-GO / REOPEN STAGE 4.**

The paper's local analytic sign-reversal mechanism survives an implementation-independent reconstruction, and the exact algebraic certificate remains valid as a narrow stationary-root/sign certificate. However, the canonical G and B3 configurations are not global public best responses on the stated strategy domain. A large profitable finite deviation moves the deviating region into a rival-public-entry / own-public-inactivity regime that the old Stage-4 production solver did not re-solve.

This is exactly the failure class that workflow v2.0 Stage 4A is designed to catch: local FOCs/SOCs and a regular-branch solver were promoted to equilibrium language without a complete all-regime finite-deviation audit.

## 2. Frozen object audited

- public choices: `x_i in [0,1]`;
- private choice: `p_T>=0` after observing public investments;
- project routes: `0,H1,H2,HT`;
- partner participation implied by `c~U[0,1]`, including boundary clipping outside the interior block;
- timing: `(x_1,x_2) -> p_T -> participation/project sorting -> surplus`;
- manuscript equilibrium language: subgame-perfect equilibrium;
- headline mathematical object: fixed-price B3 positive local public cross effect versus negative full-game G local public cross effect.

No primitive or payoff definition was modified during Stage 4A.

## 3. Independent reconstruction summary

The clean-room audit code `stage4a_v2/code/independent_adversarial_audit.py` does not import or call the Stage-4 central-regime solver or L3-1 verifier.

It reconstructs:

1. all-route project allocation as the upper envelope of primitive affine utilities on `z in [0,1]`;
2. partner masses from the primitive participation threshold with `[0,1]` clipping;
3. the participation fixed point directly from the induced demand map, with multiple starting states at audited configurations;
4. regional project surplus by direct integration of the realized upper envelope;
5. partner surplus and public cost from primitive formulas;
6. private pricing by scanning/refining the entire economically relevant positive-profit interval.

At the canonical primitives `q_T<=v+alpha=0.6`; hence for `p_T>=0.58`, `kappa_T+p_T>=0.6` and the private route cannot yield positive utility even to `z=1`. Therefore no positive-profit private optimum lies above `0.58`; the audit searches `[0,0.58]` rather than the production central-regime interval.

A separate clean-room SymPy calculation re-derives the two key beta-zero identities.

## 4. Headline theorem-certificate table

| Claim | Stage-4A result | Reason |
|---|---|---|
| First-order B3 complementarity identity | PASS, local | Independent first-order expansion gives `dM_B3/ddelta = alpha^2 d^3/(Delta_i^3 Delta_j^2)>0` |
| Full-game beta-zero cross effect | PASS, local | Independent derivation gives `M_G(0)=-3T alpha^2 kappa_L^2/[16 Delta(Delta+T)^3]<0` |
| Local strategic sign-reversal theorem | PASS, conditional/local | IFT logic is sound conditional on regular local branches and strict SOCs; no global equilibrium claim follows |
| Private-price feedback proposition | PASS as identity | It is an exact chain-rule sign condition, not a global-existence theorem |
| Exact canonical Krawczyk certificate | PASS only as isolating-box stationary-root/sign certificate | Exact local signs/curvature do not establish global public best responses |
| Canonical G public equilibrium / SPNE status | **FAIL** | Profitable finite deviation to rival-public-entry regime |
| Canonical B3 public Nash-equilibrium status | **FAIL** | Same off-regime finite deviation is profitable |
| Welfare statements premised on the canonical pair being a decentralized equilibrium | STALE / NOT CERTIFIED | Premise is false under the stated global strategy game |

See `theorem_certificates/` for the stable certificates.

## 5. Global finite-deviation / boundary / regime audit

### 5.1 Full game G

Canonical stationary configuration:

- `x_1=x_2 ~= 0.684028196361275`;
- independently re-optimized `p_T ~= 0.02274665143`;
- independently integrated `W_1 ~= 0.22301946658`.

Finite deviation by government 1:

- hold `x_2 ~= 0.684028196361275`;
- choose `x_1=0.1875`;
- re-solve downstream private price and participation across all primitive routes;
- obtain `p_T ~= 0.02518788380`;
- project masses `(n_1^F,n_2^F,n_T^F) ~= (0,0.734535,0.622508)`;
- obtain `W_1 ~= 0.25444926939`.

Gain:

`Delta W_1 ~= +0.03142980281`.

The gain is economically large relative to numerical tolerances.

### 5.2 Fixed-price B3

Canonical stationary configuration:

- `x_1=x_2 ~= 0.656020390747393`;
- matched fixed `p_T=0.02274665145`;
- `W_1 ~= 0.21091151931`.

Deviation `x_1=0.1875` gives:

- `(n_1^F,n_2^F,n_T^F) ~= (0,0.657340,0.747685)`;
- `W_1 ~= 0.23886596269`;
- gain `~= +0.02795444338`.

### 5.3 Why `x_1=0.1875` is not an arbitrary grid artifact

At that deviation, if the own hub has zero projects,

`q_1 = v + alpha(rho+x_1) = 0.26875 < kappa_L = 0.27`.

Thus `H_1` is below the outside option even for the highest project type and is genuinely inactive. The interval of own-hub inactivity contains the deviation because the threshold is `x_1<=0.19`.

Within this regime, the part of region-1 welfare that depends directly on `x_1` is

`(rho+x_1)^2/4 - gamma x_1^2/2`,

whose exact interior maximizer is

`x_1 = rho/(2 gamma - 1)=0.1875`,

with second derivative `1/2-gamma=-0.4<0`.

The deviation therefore has a transparent free-riding interpretation rather than being a numerical search accident.

## 6. Continuation audit

**FAIL.**

The previous Stage-4 routine `public_two_sided_platform_hard_kill/code/numerical_hard_kill.py` explicitly solves the smooth `0 -> H_T -> H_i` regime only. When its regime solver returns `None`, the private-profit function substitutes `-1e9`. Thus valid off-path states in other participation regimes can be treated as if they were economically impossible/unprofitable.

Workflow v2.0 requires `None`, invalid branch, nonconvergence, or regime failure to be `UNRESOLVED` unless nonexistence is separately proved. The newly found profitable rival-public-entry deviation demonstrates that the old failure semantics were economically material.

## 7. Counterexample search design and results

The audit deliberately attacked:

- public boundaries and large finite deviations over `x_i in [0,1]`;
- own-public inactivity;
- rival-public entry;
- private-price states outside the central-regime optimizer's admissible interval;
- primitive partner clipping;
- all four project route options through upper-envelope reconstruction;
- multiple participation fixed-point starting values at the audited candidate and deviation states.

A decisive counterexample was found before broader stress-grid certification was needed. Under the v2.0 kill rule, one profitable valid finite deviation is sufficient to block GO.

Permanent regression artifact: `stage4a_v2/COUNTEREXAMPLE_CANONICAL_FREE_RIDE.md`.

## 8. Benchmark-definition audit

B3 remains mathematically meaningful as a **matched-price local diagnostic**: it holds the private fee at the G stationary value while switching off repricing. The local cross-effect comparison remains useful.

What is not certified is calling the canonical B3 stationary pair a global public Nash equilibrium.

Likewise G's local stationary configuration remains useful for derivative decomposition and exact sign certification, but the phrase `subgame-perfect equilibrium` requires a globally valid first-stage public best response and valid downstream continuations after material deviations; that condition fails at the canonical vector.

No `first best` benchmark issue is implicated in this specific Stage-4A failure.

## 9. Downstream stale-state consequences

Earliest affected stage: **Stage 4**.

Because Stage 8 theory freeze relied on the canonical G/B3 equilibrium status, v2.0 rollback rules make dependent downstream states stale, including Stage 8 freeze and Stage 9--14 submission/readiness outputs, until Stage 4 is repaired/re-solved and Stage 4A passes.

The following objects are not erased and should be preserved for reuse:

- local beta-zero analytic identities;
- local sign-reversal theorem as a conditional branch theorem;
- L3-1 exact algebraic stationary-root/sign certificate;
- novelty literature work, subject to re-kill if the repaired equilibrium result changes;
- manuscript/reproducibility infrastructure.

But none of them bypasses the failed equilibrium gate.

## 10. Required repair question for Stage 4

Stage 4 must now solve the actual public game over all economically relevant project-allocation regimes and finite public deviations. In particular it must determine whether:

1. a genuine global Nash/SPNE equilibrium exists at or near the canonical primitives and still exhibits a meaningful complements-to-substitutes comparison; or
2. the research object must be reformulated as explicitly local stationary-response geometry rather than equilibrium strategic interaction.

The second route changes the equilibrium/contribution object and therefore cannot be performed silently inside Stage 4A.

## 11. Canonical verdict and routing

\[
\boxed{\textbf{NO-GO / REOPEN STAGE 4}}
\]

Reason: **profitable finite off-regime public deviation; canonical stationary pair is not a global public Nash/SPNE equilibrium.**

Do not proceed to Stage 6 under v2.0 until a repaired Stage 4 returns GO and a fresh Stage 4A returns `GO — MATHEMATICAL ADVERSARIAL CERTIFICATION PASS`.
