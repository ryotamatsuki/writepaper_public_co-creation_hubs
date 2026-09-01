# Stage 4 Contract — Minimal Model Kill Test

Date: 2026-09-01

Stage 3 preferred mechanism: **M10 — Funder-Induced Neutrality Loss / Regional Home-Bias Feedback**

Routing: **GO TO MINIMAL MODEL**, subject to the decisive prior-art/isomorphism kill test below.

## 1. Research question frozen for Stage 4

> Can a geographically bounded public sponsor's incentive to favor local matches erode an innovation intermediary's neutrality, reduce extra-regional participation, and thereby distort the realized composition of intermediation relative to a social planner in a way that is not equivalent to standard biased-intermediary or two-sided-platform theory?

Do not revert to the legacy `hub + pipeline` question.

## 2. Players

Use only the minimum set:

1. **Regional principal R** — funds/governs the intermediary and values outcomes accruing to its local constituency.
2. **One innovation intermediary I** — operates a matching/brokerage rule under fixed capacity; no profit-maximizing private intermediary benchmark at this stage.
3. **Local project/participant side L** — may be normalized as a fixed local mass if possible.
4. **Potential outside partners X** — a continuum or minimal heterogeneous-cost mass whose participation is endogenous.
5. **Social planner SP** — benchmark only, not an additional strategic player.

Do **not** add a second strategic region unless the mechanism mathematically requires it. Generic interjurisdictional competition was already killed in Stage 2.

## 3. Timing

Preferred three-stage timing:

1. Regional principal chooses or induces a local-priority / home-bias parameter `alpha ∈ [0,1]` in the intermediary's matching rule.
2. Potential outside partners observe `alpha` and choose whether to enter/participate; use a simple participation-cost cutoff if possible.
3. The intermediary matches participants according to the rule and fixed capacity; surplus is realized.

Alternative ordering is allowed only if necessary for a coherent equilibrium. Record why.

## 4. Choices

### Regional principal

Chooses `alpha`, where larger `alpha` means greater priority to local-local/local-favored matches relative to outside-partner opportunities.

### Outside partners

Choose participation `x ∈ {0,1}` or entry based on expected payoff net of participation cost. Aggregate participation `n_X(alpha)` must be derived from this decision.

### Intermediary

In the most minimal formulation, it implements the sponsor-induced priority rule rather than making a second independent optimization choice. If an intermediary best response is essential, add only one decision and explain why it cannot be folded into `alpha`.

## 5. Minimum new primitive

The only essential new primitive is:

> **Intermediary matching priority/neutrality affects an outside actor's expected chance or surplus from a match, so potential outside partners endogenously decide whether to participate.**

Do not assume directly that home bias destroys social value. The cost of bias must arise through the participant response and resulting matching-market thickness/composition.

A participation-cost distribution `c ~ F(c)` is allowed as the minimal heterogeneity required to generate an interior participation mass. Prefer a uniform or simple general CDF before introducing richer types.

## 6. Payoff discipline

Let `S_L` denote surplus accruing to local parties from local-favored matches and `S_X` surplus associated with matches involving outside partners.

The regional principal's objective may be geographically bounded, but the model **must not** obtain its headline result solely from writing

`U_R = S_L + rho S_X`, `rho < 1`,

while the planner uses `S_L + S_X`.

That would be a standard omitted-spillover result and is fatal.

Instead, `alpha` must alter `n_X(alpha)`, and `n_X(alpha)` must alter the feasible/realized matching set or matching success. The Stage 4 theorem must use this derivative/feedback essentially.

## 7. Equilibrium concept

A subgame-perfect equilibrium or backward-induction solution is sufficient:

- outsider participation best response given `alpha`;
- regional principal optimal `alpha^R` anticipating participation;
- planner optimal `alpha^SP` under the same participation technology.

No uncertainty/dynamics are required beyond participation heterogeneity unless unavoidable.

## 8. Planner benchmark

The planner chooses the same policy/matching-priority instrument `alpha` under the same technology and participation response, maximizing total surplus across local and outside actors net of participation/resource costs.

Do not give the planner a different technology.

## 9. Candidate propositions to test

### P1 — Endogenous home-bias composition wedge

Seek conditions for

`alpha^R > alpha^SP`

and therefore

`n_X(alpha^R) < n_X(alpha^SP)`.

**Kill condition:** if the proof does not use `dn_X/dalpha < 0` in an essential way, the result is generic regional weighting/spillover theory and Stage 4 fails.

### P2 — Nonlinear participation / tipping

Test whether endogenous participation creates an interior threshold or regime switch in which increasing home bias sharply reduces the viable outside pool and changes realized match composition.

A discontinuity is not required. A nontrivial threshold is sufficient if it follows from best responses rather than a fixed cost inserted solely to create the threshold.

### P3 — Local-welfare reversal

Test whether stronger local-priority pressure can eventually lower the regional principal's **own** payoff because the loss of outside participation destroys productive match opportunities.

This is the most valuable candidate result because it cannot be explained simply by the planner valuing outsiders more strongly.

### P4 — Neutrality commitment, only if P1–P3 survive

A commitment to a less home-biased matching rule may improve total welfare and possibly local welfare over a nonempty parameter region.

Do not introduce grants, subsidies, private intermediaries, or governance boards at Stage 4.

## 10. Required mathematical checks

Stage 4 must explicitly verify:

- outsider participation cutoff and monotonicity;
- existence/uniqueness or multiplicity of equilibrium as relevant;
- feasibility of matching probabilities/shares;
- regional objective FOC/SOC or global optimum characterization;
- planner FOC/SOC;
- parameter region for any strict wedge;
- whether the sign depends mechanically on an objective weight;
- boundary cases: no outside participation, full outside participation, no sponsor bias;
- comparative statics only after the main theorem survives.

Use exact algebra/SymPy where feasible and numerical diagnostics only as a supplement.

## 11. Prohibited additions

Do not add:

- two-period dynamics;
- dynamic redundancy/lock-in;
- private intermediary competition;
- multiple intermediary types;
- heterogeneous local firms beyond what is mathematically indispensable;
- political lobbying/capture as a separate game;
- digital substitution;
- endogenous absorptive capacity;
- fixed hub construction cost;
- legacy hub/pipeline variables;
- second-region strategic funding;
- KPI/multitask contracting.

If the mechanism does not work without these additions, **KILL** it.

## 12. Most dangerous prior-art comparisons

1. **de Cornière & Taylor (2019), “A model of biased intermediation.”**
2. **Zennyo (2022), “Platform Encroachment and Own-Content Bias.”**
3. canonical two-sided-market participation models summarized by **Rysman (2009)**.
4. current physical-innovation-intermediary platform framing in **Pinarello, Trabucchi & Frattini (2026)**.

The Stage 4 report must contain a model-level comparison to at least de Cornière–Taylor and Zennyo before issuing a GO verdict.

## 13. Decisive stripping/isomorphism kill test

After deriving the minimal model, remove the words:

- region;
- public;
- innovation;
- co-creation;
- hub;
- local/external.

Restate it as a generic intermediary/platform model.

If a one-to-one relabeling maps the players, bias instrument, participation response, equilibrium condition and main proposition into a standard biased-intermediary/two-sided-platform model, Stage 4 verdict must be:

**NO-GO — MECHANISM ISOMORPHIC TO PRIOR ART.**

The project survives only if regional sponsorship creates a distinct strategic term or welfare/proposition structure, not merely a new interpretation.

## 14. Stage 4 success criterion

`GO` only if all are true:

1. the minimal model is analytically coherent;
2. endogenous participation is essential to the composition wedge;
3. at least one nontrivial proposition survives;
4. the result is not driven solely by a geographically bounded welfare weight;
5. model-level comparison indicates plausible theorem-level distance from de Cornière–Taylor/Zennyo;
6. no extra features are required to rescue the mechanism.

Otherwise return `NO-GO` and recommend empirical/mixed or alternative research routes rather than repairing the legacy model.