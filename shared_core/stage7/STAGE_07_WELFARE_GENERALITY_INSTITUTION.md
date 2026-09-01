# Canonical Stage 7 — Welfare / Generality / Institutional Validation

Project: `ryotamatsuki/writepaper_public_co-creation_hubs`

Workflow: `ryotamatsuki/research-paper-workflow`

Workflow SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`

Stage 6 input branch: `shared-core/stage6-post-hardening-novelty-rekill`

Stage 6 input SHA: `34a955faa1814ac8312051b971c1f75d59aad1ee`

Stage 7 branch: `shared-core/stage7-welfare-generality-institution`

Date: 2026-09-02

---

## 1. Executive verdict

\[
\boxed{\textbf{CONDITIONAL GO}}
\]

Reason code:

\[
\boxed{\textbf{CONDITIONAL GO — REGIONAL-TO-CORE DIVERSION EVIDENCE WEAK}}
\]

The theory passes the welfare and generality tests. The only Stage-7 blocker is institutional validation of the motivating co-creation-hub application's theorem-essential diversion margin.

Do **not** proceed to Stage 7.5 until that one blocker is resolved or the co-creation-hub application is formally demoted/replaced through the workflow.

---

## 2. Welfare reconstruction

There is no separate consumer population, so consumer surplus as a separate object is **not applicable**. Welfare is aggregate realized firm/user utility net of real fixed establishment costs.

Let

\[
h=\frac{A_L-A_M-2\beta-\beta\tau}{1-2\beta},
\qquad
x=h-\tau,
\qquad
s=\frac{A_L-A_M-2\beta}{1-2\beta}.
\]

Then

\[
\boxed{SW_{00}=2C_0,}
\]

\[
\boxed{SW_{10}=SW_{01}=2C_1+\frac{h^2+x^2}{2}-F,}
\]

\[
\boxed{SW_{11}=2C_2+s^2-2F.}
\]

All three expressions were re-derived from individual utilities rather than copied from Stage 5.

---

## 3. Private versus social entry

Private/regional gross entry benefits:

\[
D_0=\beta(\tau-2h)+\frac{h^2}{2},
\]

\[
D_1=\beta(2h-2s-\tau)+\frac{s^2-x^2}{2}.
\]

Cross-region effects:

\[
X_0=\beta(\tau-2h)+\frac{x^2}{2},
\]

\[
\boxed{X_1=-\frac{\beta\tau[2-2(A_L-A_M)+\beta\tau]}{2(1-2\beta)^2}<0.}
\]

The Stage-7 prompt omitted the leading minus sign in its candidate \(X_1\) expression; the verified expression above is the correct one.

Social gross marginal entry benefits:

\[
\boxed{G_0=D_0+X_0,}
\qquad
\boxed{G_1=D_1+X_1.}
\]

The wedges are genuine cross-region welfare effects, not transfer accounting.

---

## 4. Central welfare identity

From definitions,

\[
X_1-X_0=\Gamma.
\]

Therefore

\[
\boxed{G_1-G_0=2\Gamma.}
\]

This is the most important Stage-7 welfare result.

The private strategic threshold \(A^*\) is exactly the threshold at which social marginal entry benefits switch from decreasing to increasing:

\[
A_M<A^*\Rightarrow G_1<G_0,
\]

\[
A_M>A^*\Rightarrow G_1>G_0.
\]

Hence \(A^*\) is economically meaningful beyond the private best-response game.

---

## 5. Welfare propositions

### W1 — Efficient anti-coordination: THEOREM / open-set existence

There is a nonempty open set where the asymmetric one-entry Nash equilibria coincide with the social optimum \(N^{SP}=1\).

Example with strict inequalities:

\[
\beta=0.05,\quad \tau=0.20,\quad A_L-A_M=0.60,
\]

produces

\[
D_0=0.1037654,\quad D_1=0.0838889,
\]

\[
G_0=0.1186420,\quad G_1=0.0788889.
\]

Thus \(F=0.09\) yields one-entry Nash equilibria and one-entry social optimum. Continuity gives an open set.

### W2 — Excessive proliferation under coordination: THEOREM

Define

\[
\bar G=\frac{G_0+G_1}{2}.
\]

The canonical identity is

\[
\boxed{\bar G=D_0+X_1<D_0.}
\]

In the strategic-complements region,

\[
D_0<D_1.
\]

Whenever both \((0,0)\) and \((1,1)\) are Nash equilibria,

\[
\max\{0,D_0\}<F<D_1,
\]

we necessarily have

\[
F>D_0>\bar G,
\]

so

\[
\boxed{SW_{00}>SW_{11}.}
\]

The full-entry coordination equilibrium is therefore socially excessive throughout the canonical coordination-multiplicity region.

A stronger over-entry region exists when

\[
\max\{0,\bar G\}<F<D_0:
\]

full entry is the unique Nash equilibrium while the planner chooses no peripheral entry.

### W3 — No-entry equilibrium with socially optimal full entry: FAIL / IMPOSSIBLE

No such region exists in the canonical model.

If \(\Gamma<0\), no entry requires \(F>D_0>D_1>G_1\), excluding social full entry.

If \(\Gamma>0\), no entry requires \(F>D_0>\bar G\), again excluding social full entry.

The canonical theory therefore does **not** produce a full under-entry coordination failure.

### Additional corollary — Under-provision of the first entry

In the substitutes region, the first-entry cross-region effect \(X_0\) can be positive. If \(X_0>0\), then for

\[
D_0<F<G_0
\]

no entry is decentralized equilibrium while one platform is socially optimal. Under-entry can therefore occur at the first-entry margin, but not as a no-entry-versus-full-entry coordination failure.

---

## 6. Welfare significance

The welfare content is not simply “business stealing causes excessive entry.” The residual-core mechanism determines the ordering of **both private and social marginal entry benefits**, with

\[
G_1-G_0=2(D_1-D_0).
\]

In the coordination region the planner's marginal benefits are increasing, but the negative effect of the second entry on the already-entered region implies

\[
\bar G<D_0<D_1.
\]

This produces a sharp welfare ranking: coordination can make full entry self-reinforcing even when no entry is socially preferred.

Welfare gate: **PASS**.

---

## 7. Generality assessment

Five required application families plus two stronger non-platform environments were audited.

| Environment | Result |
|---|---|
| Regional co-creation/innovation hubs | WEAK FIT |
| Local entrepreneurship/business-support hubs | WEAK / PLAUSIBLE |
| Labor-market/job-matching intermediaries | PLAUSIBLE BUT NOT EXACT |
| Local marketplaces around dominant platform | STRONG PLAUSIBLE / NEAR-EXACT |
| Research/technology-transfer intermediaries | WEAK FIT |
| Peripheral shopping malls vs downtown agglomeration | STRONG PLAUSIBLE |
| Regional airports vs dominant hub airport | STRONG PLAUSIBLE |

Retail agglomeration provides a particularly strong non-platform mapping: peripheral malls divert shoppers from downtown, while downtown attractiveness depends on variety/footfall agglomeration. Hub airports provide another: regional bypass can divert traffic from a common hub whose connectivity value depends on concentrated traffic and connection opportunities.

Generality gate: **PASS**.

---

## 8. Institutional primitive validation — co-creation hubs

| Primitive | Status |
|---|---|
| I1 Shared dominant core | SUPPORTED |
| I2 Core network value | SUPPORTED |
| I3 Regional entry diverts core participation | PLAUSIBLE BUT INDIRECT / WEAK |
| I4 Regions share the same core | SUPPORTED |
| I5 Cross-regional friction \(\tau\) | PLAUSIBLE BUT INDIRECT |
| I6 Binary establishment | SUPPORTED AS REDUCED FORM |

### I1 / I2 / I4

NEXs Tokyo explicitly uses Tokyo's concentration of information and human networks to connect startups, corporations, local governments and support institutions from across Japan. Its community, matching concierges, partners, mentors and programme structure provide strong institutional support for a common network-valued metropolitan/national core.

### I3 — binding blocker

The evidence does not establish the required diversion mechanism.

Ehime's 2026 metropolitan co-creation project explicitly seeks to strengthen the network between E:N BASE members and metropolitan companies and create mutual collaboration. Ehime's DX plan treats digital, metropolitan and local co-creation hubs as a portfolio. NEXs Tokyo's DIVE programme likewise helps regional startups use the metropolitan ecosystem.

These sources do not prove zero substitution in finite participant time/attention, but they make it unsafe to describe regional hub entry as institutionally established core depletion.

### I5 / I6

Physical local facilities, locally embedded coordinators and discrete opening decisions support \(\tau\) and binary establishment as reduced-form approximations.

Co-creation-hub application verdict:

\[
\boxed{\textbf{WEAK FIT}}
\]

Institutional gate: **one blocker remains — I3**.

---

## 9. Testable predictions

At least four externally meaningful predictions survive without theory changes:

1. the neighbor-entry effect on own establishment changes sign with core attractiveness;
2. joint peripheral-entry outcomes shift from one-entry concentration to symmetric clustering around the threshold;
3. the reversal is strongest where core value is highly installed-base/network sensitive and disappears without core-side feedback in the canonical class;
4. higher cross-peripheral friction lowers \(A^*\), because
   \[
   \frac{\partial A^*}{\partial\tau}<0;
   \]
5. social marginal entry interaction changes sign at the same \(A^*\), because
   \[
   G_1-G_0=2\Gamma.
   \]

Prediction/observability gate: **PASS**.

---

## 10. Comparative statics

Verified:

\[
\boxed{\frac{\partial A^*}{\partial A_L}=1,}
\]

\[
\boxed{\frac{\partial A^*}{\partial\tau}=-1-\frac{2\beta(1-\beta)\tau}{\Omega}<0.}
\]

For \(\beta\):

\[
\frac{\partial A^*}{\partial\beta}
=-(1-2\beta)\left[4+\frac{4\beta-16\beta^2+\tau^2}{\Omega}\right].
\]

This derivative is strictly negative under the transparent sufficient condition \(0<\beta\le1/4\). No global sign is claimed over the wider \(0<\beta<1/2\) domain.

---

## 11. Assumption classification

Essential to the canonical mechanism:

- two peripheral entrants;
- a shared residual core;
- core-side network feedback;
- endogenous user allocation;
- positive cross-peripheral friction;
- binary entry for the exact game-form result.

Tractability/minimality:

- uniform \(c\);
- linear network effect;
- symmetric regions;
- \(\beta_L=0\).

Public ownership/regional-welfare interpretation is nonessential to the core theorem.

Single-homing remains a material untested scope restriction. Stage 7 did not and may not add multi-homing.

---

## 12. Strict kill tests

### H1 — Welfare is merely standard over-entry

**SURVIVES.** Over-entry itself is standard, but the identity \(G_1-G_0=2\Gamma\) and the complete welfare ordering of the coordination region tie welfare directly to the residual-core mechanism.

### H2 — Generality is only platform relabeling

**SURVIVES.** Retail agglomeration and hub airports provide strong non-platform mappings.

### H3 — Real regional hubs do not materially thin the shared metropolitan core

**NOT RESOLVED / MATERIAL RISK.** Primary sources emphasize complementarity and network linkage; no direct diversion evidence was found.

### H4 — Strategic reversal has no organizational consequence

**SURVIVES.** The same threshold changes social marginal benefit ordering and produces systematic excessive proliferation in the coordination region.

Because H3 remains unresolved, the strict Stage-7 GO standard is not met.

---

## 13. Economic significance statement

> Decentralized entry around a network-dependent residual incumbent need not have a fixed strategic or social margin: strengthening the common core can change peripheral entry from self-limiting to self-reinforcing, and the same threshold reverses the ordering of social marginal entry benefits. This creates a sharp organizational risk in which clustered full entry can be an equilibrium even when retaining the common core without peripheral duplication yields higher aggregate welfare.

---

## 14. Referee assessment

- “No welfare content”: **MAJOR BUT ANSWERABLE**.
- “Co-creation hubs do not satisfy residual-core depletion”: **MAJOR / CURRENT BLOCKER**.
- “Regional and metropolitan hubs are complements”: **MAJOR**.
- “Only platforms”: **MAJOR BUT ANSWERABLE** via retail/airports.
- “Generality is relabeling”: **MAJOR BUT ANSWERABLE**.
- “Policy result is standard excessive entry”: **MINOR–MAJOR; concede and demote over-entry to corollary**.
- “Single-homing drives result”: **MAJOR UNRESOLVED SCOPE LIMITATION**, but not a Stage-7 authorization for model expansion.

---

## 15. Computational verification

SymPy verified all welfare identities and signs. A 200,000-draw admissible-parameter stress test found:

- \(X_1<0\) in every draw;
- \(G_1-G_0=2\Gamma\) to machine precision;
- no draw in which a no-entry Nash equilibrium coexisted with socially optimal full entry;
- every sampled coordination-multiplicity case had social optimum \(N^{SP}=0\).

No new robustness model was added.

---

## 16. Canonical verdict and routing

Canonical Stage 7 verdict:

\[
\boxed{\textbf{CONDITIONAL GO}}
\]

Reason:

\[
\boxed{\textbf{REGIONAL-TO-CORE DIVERSION EVIDENCE WEAK}}
\]

Exactly one Stage-7 blocker is authorized:

> Determine whether regional co-creation/innovation-hub participation actually substitutes at a meaningful margin for direct participation in a shared metropolitan/national ecosystem, rather than merely complementing and linking to it.

Do not proceed to Stage 7.5 while this blocker remains unresolved **if co-creation hubs are to remain the motivating/headline application**.

Permitted next action:

- targeted institutional/evidence resolution of I3 only.

Not permitted:

- changing the model;
- adding multi-homing;
- adding spillovers or congestion;
- drafting the manuscript;
- using retail/airport analogies to pretend the co-creation-hub evidence exists.

If I3 cannot be supported, the theory itself remains viable as a general IO mechanism, but the co-creation-hub application must be formally demoted or replaced before a Full-Theory Freeze Decision.
