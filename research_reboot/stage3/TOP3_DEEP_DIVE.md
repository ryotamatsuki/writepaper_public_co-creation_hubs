# Stage 3 — TOP 3 Deep Dive

Date: 2026-09-01

## TOP 1 — M10: Funder-Induced Neutrality Loss / Regional Home-Bias Feedback

### A. Minimal verbal model

A regional public principal sponsors one innovation intermediary. The intermediary connects local projects/participants with potential partners. The sponsor values successful local outcomes more strongly than surplus accruing to outside partners and can induce a matching priority/bias toward local-local opportunities. Outside partners observe the intermediary's effective priority rule before deciding whether participation is worthwhile. Because participation is costly, stronger local favoritism reduces the mass of outside participants. A smaller outside pool in turn makes external brokerage less productive, causing realized intermediation to become still more locally concentrated.

Minimum players:

1. one regional public principal;
2. one intermediary;
3. a local project/participant side;
4. a continuum or simple mass of potential outside partners with participation costs;
5. social planner benchmark.

No second strategic region is required in the minimal version. That is deliberate: the novelty test is about the sponsor–intermediary–participant feedback, not generic fiscal federalism.

### B. Mechanism diagram

`local-priority incentive of sponsor`
→ `intermediary home bias / lower outsider match probability`
→ `outside participation falls`
→ `external-side market thickness falls`
→ `realized return to external brokerage falls relative to local brokerage`
→ `equilibrium intermediation becomes more locally biased`.

The feedback must be generated through participation; it may not be inserted by assuming that external service directly has lower social/private value.

### C. Why the result is not assumed

The minimal primitive is **endogenous outsider participation in response to the intermediary's matching priority/neutrality**. The model should not contain an arbitrary reduced-form term saying that bias destroys value. Instead, an outsider of participation cost `c` enters only if its expected match surplus under priority rule `alpha` exceeds `c`. Thus a more local-biased `alpha` changes the participant pool, which changes the productivity of matching choices.

The candidate is viable only if the sponsor's best response to this participation response differs from the planner's in a non-mechanical way. If the model merely sets the sponsor's coefficient on outsider surplus below the planner's coefficient and obtains a larger `alpha`, it collapses to a generic weighted-welfare/spillover model and must be killed.

### D. Closest-paper difference

#### de Cornière & Taylor (2019)
- Existing primitive: biased advice among sellers; bias may help or harm depending on payoff congruence.
- Candidate difference to test: bias is induced by a geographically bounded public sponsor and affects whether outside innovation partners enter the intermediary's matching pool.
- Required theorem difference: sponsor-induced home bias must distort *market thickness and matching composition*, not simply re-rank existing sellers.

#### Zennyo (2022)
- Existing primitive: platform fair/biased search under first-party encroachment.
- Existing strategic feedback: biased search changes prices and consumer participation, which affects seller participation.
- Candidate difference to test: no first-party platform product or commission motive; the distortion originates in a public principal's geographically bounded objective and affects innovation-partner participation and the intermediary's service composition.
- **Danger:** this may still be mathematically isomorphic. Stage 4 must test this directly.

#### Rysman (2009) / standard two-sided platforms
- Cross-side participation is known.
- Candidate cannot claim cross-side effects as novelty.
- Required difference: regional sponsorship must create a distinct optimal-bias/composition proposition.

#### Innovation-intermediary literature
- Neutrality/trust and boundary spanning are institutionally established.
- Candidate contribution, if any, would be the strategic equilibrium/welfare consequence of sponsor-induced loss of neutrality.

### E. Candidate propositions to kill-test

**P1 — Home-bias composition wedge.** Under a nonempty parameter region, the regional sponsor chooses a more locally biased matching priority than a global planner, and the resulting outside participation is lower:

`alpha^R > alpha^SP` and `n_X(alpha^R) < n_X(alpha^SP)`.

This proposition is acceptable only if the wedge depends on the participation response and not solely on a lower sponsor weight on outsider surplus.

**P2 — Participation tipping / corner.** If participation costs are heterogeneous and expected match probability falls sufficiently with home bias, there may be a threshold `alpha_bar` beyond which the outside side becomes too thin for external matching to be sustained. The existence of a discontinuity is not required; an interior nonlinear threshold is sufficient.

**P3 — Non-monotone local objective.** Stronger local-priority pressure can initially redirect matches toward local actors but eventually reduce even the sponsor's local payoff because lost outsider participation destroys valuable complementary match opportunities. This is the most attractive theorem shape because the cost of bias is generated endogenously through the matching pool.

**P4 — Neutrality commitment (optional).** If P1–P3 survive, a commitment to a less biased matching rule may improve both total welfare and, over some region, the sponsor's own payoff. Do not add a subsidy or governance instrument unless the basic mechanism works first.

### F. Expected failure modes

Kill M10 if any of the following occurs:

1. `alpha^R > alpha^SP` follows solely because the sponsor mechanically ignores outsider surplus.
2. outsider participation is unnecessary for the result;
3. stripping labels makes the model a direct special case of de Cornière & Taylor or Zennyo;
4. the only result is monotonic `more local weight → more local bias`;
5. the matching technology must contain an imposed complementarity to generate the feedback;
6. the model needs dynamics, heterogeneous local firms, private intermediaries, or political capture to work.

### G. Stage 4 complexity

Expected minimal system: 2–4 first-order/threshold equations: outsider participation cutoff, matching/output expression, regional sponsor choice of `alpha`, planner choice of `alpha`. Tractability: **HIGH–MEDIUM**.

---

## TOP 2 — M8: Disclosure-Mediated External Search with Endogenous Participation

### A. Minimal verbal model

An intermediary chooses how much identifying/project information to reveal to potential outside partners. More disclosure can improve search and matching but may impose appropriation/confidentiality costs on local innovators. Local participants and outside solvers decide whether to participate based on the disclosure regime.

### B. Mechanism diagram

`disclosure/anonymization rule`
→ `local and outside participation`
→ `matching thickness/quality`
→ `value of disclosure service`
→ `disclosure rule`.

### C. Why the result might be endogenous

Participation choices respond to disclosure rather than an exogenous `pipeline benefit`. The service composition could tilt toward local/private search when disclosure costs are high and toward broad external search when anonymity or protection is effective.

### D. Closest-paper difference and threat

Kawakami (2024) already gives a competitive matching model with professional disclosure services, endogenous disclosure, monopoly provision, a minimum-disclosure option and welfare analysis. Boudreau & Lakhani provide innovation-specific disclosure/incentive theory. Howells & Thomas describe intermediary search, selective revealing and confidentiality roles.

Therefore the generic M8 mechanism is **very close to prior art**. It remains TOP 3 only because regional/public objectives could potentially add a distinct composition problem, but that would risk feature accumulation.

### E. Candidate propositions

- disclosure can shift the local/external participant mix non-monotonically;
- a public disclosure minimum can improve participation composition under some parameters;
- excessive disclosure can reduce participation by high-value local projects.

### F. Failure modes

Kill if public/regional sponsorship must be added as a second independent mechanism to obtain novelty. Kill if it is a direct specialization of Kawakami's disclosure-service model.

### G. Stage 4 complexity

MEDIUM. Information structure and selection create more moving parts than M10.

**Status: SECONDARY.**

---

## TOP 3 — M3: Service-Mix-Induced Participant Sorting / Market-Thickness Feedback

### A. Minimal verbal model

A multitask intermediary allocates effort between local and extra-regional matching. Different service mixes attract different masses/types of participants. Entry then changes the return to the two services, creating a feedback between service composition and the composition of the matching pool.

### B. Mechanism diagram

`service mix L/X`
→ `local/outside participation`
→ `relative market thickness`
→ `match success on L/X`
→ `service mix L/X`.

### C. Why the result is not necessarily assumed

The feedback can arise from participation cutoffs and matching probabilities rather than from an imposed complementarity between `L` and `X`.

### D. Closest-paper difference and threat

The mechanism sits directly inside canonical two-sided-platform logic. Rysman makes cross-side participation standard; 2026 work also explicitly describes physical innovation intermediaries as two-sided platforms. Carayol & Sterzi already make intermediary use endogenous. Thus M3 is not a defensible contribution unless the public/regional institution creates an additional *single* primitive that changes the theorem.

### E. Candidate propositions

- multiple service-composition equilibria or tipping may arise through participation feedback;
- decentralized service composition may differ from planner allocation;
- a small change in one side's participation cost can trigger a discrete or strong shift in the service portfolio.

### F. Failure modes

Kill if it is only a two-sided platform with `local` and `external` replacing `side A` and `side B`. Do not rescue it by adding a public principal, dynamics and heterogeneity simultaneously.

### G. Stage 4 complexity

HIGH tractability, but low closest-paper survival relative to M10.

**Status: SECONDARY.**

---

# Preferred candidate

**M10 — Funder-Induced Neutrality Loss / Regional Home-Bias Feedback.**

This choice is not based on score alone. M10 has the clearest chance to connect a real institutional feature of public innovation intermediaries—dependence on a geographically bounded sponsor—with a productive attribute of intermediation—neutrality/trust that sustains participation by outside actors.

The Stage 4 burden is deliberately severe: show that this interaction produces a theorem that does not reduce to standard biased intermediation, platform self-preferencing, or a regional planner putting a lower weight on outsiders. If that cannot be shown in the minimal model, return `NO-GO` rather than add features.