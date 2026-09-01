# Stage 4 — Minimal Model Construction and Isomorphism Kill Test

## Public Co-Creation Hubs — Autonomous Research Reboot

Date: 2026-09-01

Stage 4 base SHA: `a4bff663b34898fd61a7eb67d4ac989515b5aca5`

Workflow SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`

Branch: `research-reboot/stage4-minimal-neutrality-model`

## 1. Executive verdict

**Canonical verdict: NO-GO.**

Specific classification:

**NO-GO — RESULT TOO WEAK.**

Routing under the canonical workflow: **terminate this theory branch.**

The Stage 3 preferred mechanism does not collapse mathematically. A minimal model produces:

- endogenous outsider participation falling with home bias;
- a unique interior regional own-welfare-maximizing bias;
- a planner choosing less bias;
- a nonempty open parameter region;
- a local-welfare reversal in which excessive favoritism harms the sponsor itself;
- disappearance of that local-welfare reversal when participation is fixed.

However, the decisive novelty gate fails. Once institutional labels are removed, the model is a generic biased-principal/intermediary participation trade-off. Public/regional sponsorship does not generate an additional strategic term, state variable, information problem, or best response. Replacing the regional sponsor with a private biased platform/principal leaves all equations and propositions unchanged.

The model is not literally one-to-one isomorphic to de Cornière & Taylor (2019) or Zennyo (2022), so `NO-GO — MECHANISM ISOMORPHIC TO PRIOR ART` would be too strong. But the surviving theorem is not sufficiently distinct from the established biased-intermediation/platform-participation literature to justify formal development.

## 2. Binding Stage 3 contract

Stage 3 authorized only M10:

> Funder-Induced Neutrality Loss / Regional Home-Bias Feedback.

The model was frozen to one sponsor/principal, one intermediary, a fixed local/favored side, endogenous outside participation, one bias parameter, and a planner benchmark.

No dynamics, second strategic region, private intermediary, KPI contract, absorptive capacity, legacy hub/pipeline technology, network formation, or other rescue feature was added.

## 3. Exact minimal mechanism

The home-priority parameter is

`alpha in [0,1]`.

A unit mass of potential outside participants has

`c ~ Uniform[0,1]`.

Given outside participation `n`, entrant gross expected access value is

`g(alpha,n)=v(1-alpha)/(1+n)`.

The denominator represents congestion under fixed intermediary attention capacity.

Outside participation solves

`n=g(alpha,n)`,

hence

`n(alpha)=[sqrt(1+4v(1-alpha))-1]/2`.

For `0<v<2`, the CDF upper bound never binds.

## 4. Players and timing

1. Sponsor chooses `alpha`.
2. Outside participants observe `alpha` and enter using a participation-cost cutoff.
3. The intermediary implements the priority rule and the participant pool determines brokerage/search value.
4. Planner benchmark uses the same technology.

The intermediary has no second optimization problem, in compliance with the Stage 3 minimality contract.

## 5. Participation problem

Define

`D=1+4v(1-alpha)`.

Then

`n(alpha)=[sqrt(D)-1]/2`,

`n'(alpha)=-v/sqrt(D)<0`,

and

`n''(alpha)=-2v^2/D^(3/2)<0`.

Thus more home bias reduces outside participation, and the marginal participation loss becomes larger as bias rises.

## 6. Regional objective

Regional own payoff is

`U_R(alpha)=A alpha+B n(alpha)`.

`A>0` is the marginal local value of home-priority intermediation.

`B>0` is the local value of breadth/quality in the outside opportunity pool assembled by the intermediary.

The key discipline is that `U_R` is linear in `alpha` and `n`; no ad hoc quadratic local-welfare loss is inserted.

Regional derivatives are

`U_R'=A-Bv/sqrt(D)`

and

`U_R''=-2Bv^2/D^(3/2)<0`.

Therefore the regional problem is globally strictly concave.

## 7. Regional equilibrium

The interior regional solution is

`alpha_R=1+1/(4v)-B^2 v/(4A^2)`.

Equivalent participation is

`n_R=(vB-A)/(2A)`.

The regional optimum is interior iff

`Bv/sqrt(1+4v)<A<Bv`.

## 8. Planner objective

At participation equilibrium, the marginal entrant's gross value equals `n`. Aggregate outside net surplus is

`S_X=n^2/2`.

Therefore

`W(alpha)=A alpha+B n(alpha)+n(alpha)^2/2`.

The same matching/participation technology is used.

Derivatives are

`W'=A-v/2-v(B-1/2)/sqrt(D)`

and

`W''=v^2(1-2B)/D^(3/2)`.

For `B>1/2`, planner welfare is strictly concave.

## 9. Planner optimum

The interior planner solution is

`alpha_SP=1+1/(4v)-v(B-1/2)^2/[4(A-v/2)^2]`.

Equivalent participation is

`n_SP=(vB-A)/(2A-v)`.

A joint open interior region is

`Omega={0<v<2, B>1/2, v/2+v(B-1/2)/sqrt(1+4v)<A<vB}`.

## 10. Composition wedge

At the regional optimum,

`U_R'(alpha_R)=0`.

Since

`W'=U_R'+n n'`,

we obtain

`W'(alpha_R)=n(alpha_R)n'(alpha_R)<0`.

Strict concavity of `W` then gives

`alpha_SP<alpha_R`

and therefore

`n_SP>n_R`.

So the regional sponsor chooses more home bias and supports a thinner outside pool than the planner on `Omega`.

This is algebraically clean but is not by itself a distinct contribution.

## 11. Participation threshold

The minimal model does **not** generate an interior tipping threshold.

`n(alpha)>0` for every `alpha<1`, with

`n(1)=0`.

Creating an interior threshold would require a new fixed cost, minimum cost, coordination primitive, or discrete-entry structure. None is authorized.

Candidate P2 is therefore killed.

## 12. Local-welfare reversal

This is the strongest mathematical result.

On the regional interior region:

`U_R'(0)>0`,

`U_R'(1)<0`,

and

`U_R''<0`.

Hence there is a unique interior regional optimum, and for

`alpha>alpha_R`,

`U_R'(alpha)<0`.

Thus excessive home favoritism reduces even regional own welfare because it drives away participants whose presence is locally valuable.

This is not created by the planner's broader welfare weight.

## 13. Fixed-participation counterfactual

Set `n=nbar` exogenously.

Regional payoff becomes

`U_R^F=A alpha+B nbar`.

Therefore

`dU_R^F/dalpha=A>0`

and

`alpha_R^F=1`.

The local-welfare reversal disappears completely.

This establishes that P3 is genuinely participation-essential inside the minimal model.

However, a planner/regional wedge can still survive under fixed participation because outsider direct payoff falls with bias. Thus P1 alone can be generated by a generic welfare-weight difference and is not sufficient.

## 14. Symbolic verification

`stage4_symbolic_verification.py` verifies from primitives:

- participation fixed point;
- `n'` and `n''`;
- regional FOC/SOC;
- regional closed form;
- planner FOC/SOC;
- planner closed form;
- `W'=U_R'+n n'`;
- participation identities at both candidate optima;
- fixed-participation derivatives.

Status: **PASS**.

## 15. Numerical parameter-region validation

Benchmark:

`v=1, B=1, A=0.8`.

Then

`alpha_R=0.859375`,

`alpha_SP=0.555555...`,

`n_R=0.125`,

`n_SP=0.333333...`.

Regional own payoff satisfies

`U_R(alpha_R)=0.8125 > U_R(1)=0.8`.

The numerical audit also checks 25 admissible points spanning five values of `v` and five values of `B`; all satisfy the expected interior ordering and welfare inequalities.

This is a sanity check only; analytic proofs are available.

## 16. Proposition audit

### P1 — Home-bias composition wedge

**TRUE, but KILL AS CONTRIBUTION.**

The exact ordering exists, but a regional/planner wedge can survive when participation is fixed because of standard welfare-weight differences.

### P2 — Participation tipping

**FALSE in the minimal model.**

No interior tipping point exists.

### P3 — Local-welfare reversal

**TRUE AND PARTICIPATION-ESSENTIAL, but NOVELTY FAILS.**

This is the strongest result mathematically.

### P4 — Neutrality commitment

**NOT PURSUED.**

A distinct commitment result would require another mechanism such as agency, time inconsistency, or governance design.

## 17. Prior-art isomorphism / stripping test

After removing all institutional labels, the model is:

> A principal chooses a bias parameter. Bias raises direct value from a favored side but reduces participation on another side. The principal values that second-side pool, so excessive bias can lower its own payoff. A planner values participant surplus additionally and chooses less bias.

This is not literally identical to the exact models of de Cornière & Taylor or Zennyo:

- de Cornière & Taylor emphasize strategic seller responses to intermediary advice;
- Zennyo includes vertical integration, commissions, prices, and two-sided participation;
- Stage 4 uses one participation side and congestion.

But these differences do not deliver a new mechanism.

The economically active Stage 4 loop is the standard class:

`bias -> participation response -> intermediary/principal payoff -> optimal self-limiting bias`.

Current literature makes the problem more severe, not less. Platform-bias theory now includes equilibrium welfare measurement, self-preferencing with multiple instruments, consumer search, and nonlinear biased matchmaking.

## 18. Innovation-intermediary institutional interpretation

The application is credible.

Existing innovation-intermediary research explicitly identifies:

- neutrality and impartiality;
- trust and credible positioning;
- governance/funding tensions;
- publicly funded neutral matchmaking.

Accordingly, the model has an empirical/institutional bridge.

But this does not rescue theorem-level novelty: neutrality is already an established institutional property, while the formal bias-participation discipline is established platform/intermediation economics.

## 19. Strongest referee attack

> “Replace the regional public sponsor by any private principal that benefits from favoritism and also values participation. Nothing in the model changes. The local-welfare reversal is a standard demand/participation discipline on bias, not a new theory of regional innovation intermediation.”

Assessment: **FATAL.**

## 20. Why the result is too weak rather than algebraically false

The Stage 4 mechanism passes the internal mathematical tests:

- exact solution;
- unique interior region;
- participation essentiality for P3;
- local-welfare reversal;
- meaningful boundaries.

The failure is external validity as a contribution.

The regional/public sponsorship interpretation does not create a distinct strategic object. Making sponsorship economically essential would require adding governance, agency, common-agency funding, political accountability, or another mechanism. Those changes are prohibited and, in several cases, already faced strong prior art in Stage 3.

Therefore the model must not be repaired at this stage.

## 21. Canonical Stage 4 verdict

**NO-GO**

Specific diagnosis:

**NO-GO — RESULT TOO WEAK.**

This is preferable to `NO-GO — MECHANISM ISOMORPHIC TO PRIOR ART` because there is no literal one-to-one mapping to either de Cornière–Taylor or Zennyo. The problem is structural genericity rather than exact duplication.

## 22. Next-stage contract

Under the canonical workflow, `NO-GO` terminates this branch.

Do **not** proceed to Stage 5 or Stage 6 for M10.

Do **not** add governance, a second region, private intermediaries, dynamics, KPI incentives, absorptive capacity, or network formation to rescue M10.

Permissible future routes require a genuinely distinct pivot:

1. return to Stage 3 with a different mechanism not already killed, only if there is a specific reason to reopen mechanism search; or
2. return to Stage 0 with a different research question; or
3. pivot the co-creation-hub project toward empirical/mixed institutional analysis, where neutrality, funding structure, local/outside participation composition, and governance can be measured rather than claimed as a new general theorem.

The legacy paper remains frozen and is not approved for revision/submission.