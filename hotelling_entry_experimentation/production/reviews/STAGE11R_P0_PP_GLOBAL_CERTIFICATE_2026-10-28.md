# Stage 11R / targeted Stage 4A rollback — exact P/P fee equilibrium and four-actor entry certification

**Date:** 2026-10-28 JST. **Frozen theory:** \`STAGE33-2026-10-28-THEORY-FREEZE-v1\`. **Workflow basis:** \`research-paper-workflow/THEORY_PAPER_RESEARCH_PIPELINE.md\`, Stage11 regression to Stage4/4A, and cross-stage §3.2, §3.11. **Decision:** **P0 BASELINE MATHEMATICAL REPAIR PASS / STAGE11 OVERALL HOLD**. The word PASS is **only** for the three exact numerical baseline environments listed below. No parameter-global generality, new theory, physical-hub superiority, external-human verification, or actual journal-submission permission.

## 1. Discovery / provenance / change control

Earlier \`production/verification/stage33_reproduce.py\` used sampled fees and damped best-response iteration for off-path \`(P,P)\` continuation, and explicitly disclaimed global pure-Nash verification. The separate \`production/reviews/stage11_independent_entry_audit.py\` reproduced *entry flags from already exported CSV*, but did not construct fee continuations. Stage11 thus declared an honest **MAJOR certification regression**, routed to earliest affected mathematical stage rather than patched in manuscript prose.

This audit derives a **new certificate of existing, unchanged Stage8 baseline**, without replacing its primitives. The new \`production/reviews/stage11r_pp_global_certificate.py\` does **not** import original fee solvers. It uses \`fractions.Fraction\` for exact rational arithmetic, rational endpoint interval arithmetic, Poincaré–Miranda face signs, and independently reconstructed continuation/entry payoffs. Its root-center decimals originate from a numerical pilot, but all *claimed strict inequalities* are verified over rational enclosing intervals, not presumed from floating roots. The pilot search is neither the theorem nor its mathematical proof. Verification actors: **AI** (mathematical argument), **COMPUTATION** (rational checks), **EXTERNAL HUMAN: none**, **AUTHOR: final approval pending**.

## 2. Mathematical lemma: global P/P price Nash in the specified baseline

Assume two \`P\` providers, region \`i=1,2\`, population \`N_i>0\`, common \`Delta_i=1\`, uniform founder valuation equal to technical ability, and fixed frozen nonnegative coefficients \`u_i>0,a_i>=0,b_i>=0\`. At posted nonnegative fees \`f_i\`, let \`x_i\` be the participation cutoff (not the founder ability or firm-design variable from another paper). The coupled cutoff fixed point is:

\`\`\`text
x_i = clip[ (f_i-u_i-a_i(1-x_j)) / D_i(x_j) , 0, 1 ]
D_i(y)=1+(b_i/2)(1-y^2),   k_i=a_i+b_i.
\`\`\`

At interior cuts, the derivative of provider i's cutoff response w.r.t. the opposing cutoff is \`(a_i+b_i*x_i*x_j)/D_i <= k_i\`; at clip branches the derivative is zero and continuity preserves the same global Lipschitz bound. Since the exact baseline checks \`max(k_A,k_B)<1\`, the simultaneous cutoff map is a sup-norm contraction for **every** nonnegative fee vector. **Thus participation exists and is unique for every fee vector**; it is not just an artefact of damped iteration.

Fix opponent fee \`f_j\`. For any potential own cutoff \`x in [0,1]\`, first calculate the opponent cutoff uniquely and explicitly as \`y(x) = clip[(f_j-u_j-a_j(1-x))/D_j(x),0,1]\`. The own inverse fee is

\`\`\`text
f_i(x;f_j)=u_i+a_i(1-y(x))+D_i(y(x))*x.
y'(x)= (a_j+b_j*x*y)/D_j(x)  if 0<y<1; otherwise 0.
f_i'(x)=D_i(y) - (a_i+b_i*x*y)*y'(x).
\`\`\`

Every derivative (including one-sided kink derivatives) satisfies \`1-k_i*k_j <= f_i'(x) <= 1+b_i/2\`; therefore the inverse fee is strictly increasing since \`k_A*k_B<1\`. Fees below \`f_i(0;f_j)\` yield full participation and revenue increasing in fee, so they cannot be an optimum. Fees above \`f_i(1;f_j)\` generate zero revenue, and are also never optimal. Hence all globally optimal fees are represented by maximising \`R_i(x)=N_i*f_i(x;f_j)*(1-x)\` on \`x in [0,1]\`.

Using \`u_i+x <= f_i(x;f_j) <= u_i+a_i+(1+b_i/2)x\` and the derivative bounds, each global best response cutoff is strictly inside \`[L_i,U_i]\`, where

\`\`\`text
z=k_A*k_B
L_i=(1-z-u_i-a_i)/(2+b_i/2-z)
U_i=(1+b_i/2-u_i)/(2+b_i/2).
\`\`\`

Indeed \`R_i'(x)>0\` for \`x<L_i\` and \`R_i'(x)<0\` for \`x>U_i\`; one-sided derivatives cover kinks. Every best-response fee is thus in \`K_i=[u_i+L_i, u_i+a_i+(1+b_i/2)*U_i]\`. **The exact rational audit checks** \`0<L_i<U_i<1\` and \`u_i+a_i<min K_i<=max K_i<u_i+1\` separately for both providers in all three baselines. It follows that for any opponent fee in \`K_j\`, \`0<y(x)<1\` for *all* \`x in [0,1]\`: no opponent mass-entry or mass-exit kink is omitted on the entire global own-deviation path.

For opponent fees in \`K_j\`, \`y''=b_j(y+2*x*y')/D_j>=0\` and

\`\`\`text
f_i''(x)= -b_i*(2*y+x*y')*y' - (a_i+b_i*x*y)*y'' <= 0
R_i''(x)=N_i*((1-x)*f_i''(x)-2*f_i'(x)) < 0.
\`\`\`

Thus every in-box opposing fee gives a unique **global** own-cutoff best response; no numerical grid, local FOC-only check, or silent exclusion of full-coverage/off-path fee deviations is needed. Set \`H_i=b_i+b_i*k_j+b_j*k_i+k_i\`. Since \`|partial f_i / partial f_j|<=k_i\`, \`|partial f_i'/partial f_j|<=b_i+b_i*k_j+b_j*k_i\` and \`partial R_i'/partial x<=-2N_i(1-z)\`, the posted-fee best-response map satisfies

\`\`\`text
|d BR_i(f_j) / d f_j| <=
  k_i + (1+b_i/2)*H_i/[2*(1-z)] =: B_i.
\`\`\`

The code verifies \`B_i<1/100\` by **exact fractions** for each provider and each baseline. The BR map therefore sends the closed rectangle \`K_A x K_B\` continuously into itself and contracts in the sup norm. Banach's theorem gives **existence and uniqueness of a pure P/P fee Nash equilibrium** on this rectangle; every nonnegative-fee pure Nash must be within the same rectangle because every global BR is in \`K_i\`. Consequently this is **global pure Nash existence and uniqueness among all nonnegative posted prices** for each audited baseline, not a claim about arbitrary parameter vectors.

This is an analytic sufficient-condition theorem with exact baseline checks. It is **not** a new economics research contribution or a fully formalised Lean proof.

## 3. Rational location certificate of the unique fee equilibrium

Let \`D_i=1+b_i*(1-x_j^2)/2\`, \`C_i=a_i+b_i*x_i*x_j\`, \`f_i=u_i+a_i*(1-x_j)+D_i*x_i\`, and

\`\`\`text
F_i(x_A,x_B)= (1-x_i)*(D_i-C_i*C_j/D_j)-f_i.
\`\`\`

At the unique price Nash, both \`F_i=0\`. Independently of floating-point optimisers, the checker evaluates the entire **rational** interval extension of \`F_A\` at the lower and upper A-faces of the box, and similarly \`F_B\` on B-faces, each time allowing the other cutoff to vary across its full box interval. The four signs are strictly \`(+,-,+,-)\`. The continuous Poincaré–Miranda theorem gives an interior simultaneous stationary point; global concavity and BR contraction proved above establish that it is the **unique global** P/P fee Nash.

In all three rows the cutoff box radius is \`1/1,000,000,000\` about the following exact decimal rational centres (not rounded point-equality assertions):

| lambda | technology | x_A centre | x_B centre | fee_A approximately | fee_B approximately | max B_i (upper bound) |
|---:|---|---:|---:|---:|---:|---:|
| 0 | ability sorting | 0.380848479444378 | 0.359890064348462 | 0.619151416551 | 0.640109828126 | 0.000709887 |
| 0.5 | ability sorting | 0.381077381587935 | 0.360088394043391 | 0.620114272513 | 0.640818991665 | 0.009589462 |
| 0.5 | pooled mean-success | 0.380493572505174 | 0.359632374104520 | 0.619505267789 | 0.640366427138 | 0.002369802 |

All reported fee approximations above lie inside *strict rational inclusion enclosures* stored in the generated JSON. All 24 exact face-sign checks across lower/upper faces and region/scenario were run; there are **12** actual faces: 4 per scenario. Do not misread the loop count as 24 faces.

## 4. First-stage full revalidation and smallest slack

With the P/P fee equilibrium now uniquely established at all off-path prices in the three baselines, the checker independently reconstructs all **16 entry profiles** per scenario (no exported 48-row CSV), all **four one-player binary government/provider deviations** per profile, and all **15 competing restricted-planner alternatives**. Every alleged NE condition has an interval with strictly positive lower bound. For every other configuration, at least one strictly profitable deviation has been proven by an interval with strictly negative upper bound. For each planner competitor, the \`(B,B)\` welfare gap has a strictly positive rational lower endpoint.

| lambda | technology | Unique first-stage pure NE | Unique restricted planner | Minimum selected NE slack lower bound | Minimum restricted-planner gap lower bound |
|---:|---|---|---|---:|---:|
| 0 | ability sorting | (B,P) | (B,B) | 0.0000749289928794 | 0.000338879255049 |
| 0.5 | ability sorting | (B,B) | (B,B) | **0.0000316701024488** | 0.000799326991258 |
| 0.5 | pooled mean-success | (B,P) | (B,B) | 0.000456770766838 | 0.000936816106723 |

The \`0.0000316701\` near-tie survives the exact rational interval proof. Its fragility in parameter space, the narrow public cost band \`(0.040043229233,0.040531670102)\`, and the adverse Beta/Bertrand scope warnings remain genuine and must NOT be erased or called robust/generic.

Two distinct statements are now established for the three frozen baseline environments only: **(i)** all material P/P continuations exist uniquely and obey global tariff best responses, and **(ii)** the induced four-player entry equilibria and restricted planner rankings pass strict rational interval audits. Original manuscript only had numeric-gridded version of (i). Generic entry equilibrium existence, Beta(2,2) global fee uniqueness, all-parameter government-ownership thresholds, cross-conduct portability and an unconditional public physical hub case are not shown.

## 5. Reproduction and evidence labels

Run in an independent directory:

\`\`\`bash
cd hotelling_entry_experimentation/production/reviews
python3 stage11r_pp_global_certificate.py
# writes stage11r_pp_global_certificate.json
\`\`\`

This uses only Python 3 standard library. Every underlying parameter is an exact \`Fraction\`. Centers are exact decimal rationals. \`I\` does endpoint-fraction inclusion arithmetic, so strict comparisons do not depend on the printed approximate values. The JSON logs 3 scenarios × (4 face signs + 16-profile response classification + planner comparison), and each scenario's exact-checked price-contraction box/bounds. The output JSON itself prints decimal endpoints for readability; compare decisions in the exact checker, not rounded displays. If GitHub Actions has not yet run on this PR head, **do not claim a green remote CI**; local independent Python 3.13 execution produced all PASS.

No frozen Stage8 source, certificate, Stage10 final manuscript, Stage11 original reviewer report or Stage12 journal ranking was edited. Existing Lemma/Stage4A claims were not silently enlarged. This new certified auxiliary fact should be incorporated as a narrowly scoped P/P global equilibrium proof into the manuscript appendix only after an independent referee-verifiability review.

## 6. Stage gate / next action

- **Targeted Stage4A correctness/globality regression:** **PASS for all three baseline scenarios**, from the analytic lemma plus exact rational certificate, not only a dense finite grid.
- **Stage11 mathematical continuation defect:** **REMEDIED in the baseline scope**; manuscript appendix update and independent second review pending.
- **Stage7.5A quantifier portability:** **unchanged** — no transfer to Beta, arbitrary fee coefficients or different product competition.
- **Stage11 scientific novelty/value:** **HOLD** — P1 and P3 remain known specializations; P4 is a narrow but now mathematically validated numerical existence and can be submitted only with proper scope.
- **Stage12 editorial journal target:** provisional JICT fit unchanged. Do **not** let journal preference trigger model changes.
- **Stage13:** not formally closed. Next bounded work item after independent scrutiny of this certificate: integrate the global P/P lemma, refocus Abstract/Introduction around the integrated economic mechanism, and perform direct recent-prior comparison and bounded empirical predictions. Do not reopen Stage8 solely to improve journal prestige.

**Caution:** The Python script plus mathematical derivation is not a proof-assistant formalisation, not five external human referees, and not author personal proof verification. A genuine new mathematical counterexample should reopen Stage4A rather than be hidden under an unearned PASS.
