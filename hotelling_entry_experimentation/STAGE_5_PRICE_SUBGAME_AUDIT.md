# Stage 5 (2026-10-27): Price continuation multiplicity and SPNE gate

## Verdict
**BLOCKER IDENTIFIED.** The Stage 4 interior FOC/SOC candidate is valid as a conditional local calculation, but its global SPNE and unconditional private-entry/municipal-establishment claims are NOT certified. At off-path locations where the private service dominates socially, the public welfare-maximizer's price is payoff-irrelevant over a continuum, producing multiple pure price Nash equilibria with sharply different private profits.

This is an audit of the inherited model, not a novel publication theorem.

## Restricted two-hub price subgame
Mass M of types x uniform on [0,1], public characteristic a, private characteristic b, quadratic mismatch t>0, private direct service advantage delta>0:
U_G(x)=v-p_G-t(x-a)^2;
U_P(x)=v+delta-p_P-t(x-b)^2.
Nonparticipation utility=0, p_G,p_P >=0. Public maximizes total gross surplus (consumer surplus plus firm profits, no distortionary public funds); private maximizes p_P times demand. All consumers are fully covered in the displayed examples.

For a<b write D=b-a>0. Efficient consumer split (when interior):
z_eff=(a+b)/2-delta/[2tD].
If 0<z_eff<1 and no outsiders, public optimal pricing implies p_G=p_P and private price FOC yields p=2tD*(1-z_eff).
These are the Stage 4 expressions. They do NOT characterize all price regimes.

## Lemma: continuum of price equilibria in all-private regime
Suppose
Delta := delta - tD(a+b) > 0,
so that service P (ignoring the transfer price) is socially superior to G for every x in (0,1].
Assume nonnegative prices and select q satisfying
  max{0, 2tD-Delta} <= q < v-t*max{a^2,(1-a)^2}.
Then the following is a Nash equilibrium of the price subgame:
  p_G=q, p_P=q+Delta.
At these prices, P serves all users except possibly the zero-mass type x=0; everyone is covered by the reservation utility bound.

Proof:
  U_P(x)-U_G(x)=2tDx >=0.
Public total welfare is globally maximized because full participation is achieved and every type is assigned to the socially superior provider P; any public price that accomplishes this is a best response.
Given q, for p_P<=q+Delta the private hub serves the entire market and prefers the highest price. For p_P>q+Delta, private demand is 1-(p_P-q-Delta)/(2tD), truncated at zero, while G at q can cover everyone. Private profit immediately to the right of the full-market kink is nonincreasing iff q+Delta>=2tD, and is strictly decreasing thereafter. Thus p_P=q+Delta is a best response.
As q ranges over a nondegenerate interval, private equilibrium profit M*(q+Delta) varies. QED.

The case b<a is symmetric with Delta_left=delta-t(a-b)(2-a-b). Near exact collocation, a tie-breaking convention is necessary; the strict b>a construction above does NOT rely on positive-mass ties.

## Verified cases (t=1,v=2,M=1, original Stage 4 calibration)
Case n=1: delta=.09, a=.16, original b=.66; conditional interior profit .4624. Deviate to b=.17; D=.01, Delta=.0867, q=1, private price1.0867, private profit1.0867. Gain=.6243.
Case n=2: delta=.15, a=.10, original b=.60; conditional interior profit .6400. Deviate to b=.11; D=.01, Delta=.1479, q=1, private price1.1479, private profit1.1479. Gain=.5079.
Price examples satisfy strict full coverage and all best-response inequalities. A 30,001-point price search also verified the private price is profit-maximizing over [0,3]. These off-path price equilibria demonstrate failure of "deviation robustness to every price equilibrium continuation" — not failure of all possible SPNE profiles.

## A LOW-price continuation can still support the candidate in a restricted search
For all-private-dominance locations, choose the minimum viable public price:
  q_low=max(0,2tD-Delta).
The associated private profit is max(Delta,2tD).
For interior-split locations, use p_G=p_P=2tD*(private share). For all-public dominance, private profit zero.

With this explicit continuation-selection convention and the same v=2, t=1 benchmark, dense 10,001-point searches of private position b across [0,1] find the original b=.66 (delta=.09,a=.16) and b=.60 (delta=.15,a=.10) maximize PRIVATE profits globally over the tested branch schedule. Independent 100,001-point grid checks of the socially efficient public allocation show a=.16 and a=.10 maximize public welfare at their original opponent positions.

IMPORTANT: THESE ARE GRIDS, NOT AN EQUILIBRIUM EXISTENCE PROOF. For a complete pure-strategy SPNE one must specify Nash prices at every location profile, prove continuation equilibria for each feasible location and participation regime (including outside option), handle location crossovers and collocation, then rule out all location and establishment deviations. The model may require a payoff-relevant refinement or explicit government price commitment for determinate comparative statics.

## Consequences for prior stages
- Previous Stage 4 "private profit = 4 x marginal social benefit" is an algebraic identity at the identified INTERIOR candidate only. It is NOT a full equilibrium theorem.
- Previous private-entry thresholds and municipal build/no-build schedule are CONDITIONAL on an equilibrium-selection convention and interior applicability.
- The Stage 1 position numerics lack an archived optimizer and remain unverified.
- Do NOT claim that the full game has no SPNE: price-equilibrium multiplicity by itself does not imply SPNE nonexistence.
- Do NOT claim unique entrant profits, or unique build thresholds, without specifying off-path continuation and selection.

## Research decision
Next minimal step: finish price subgames piecewise with nonparticipation (especially for post-deviation prices), establish a canonical equilibrium-selection or economically microfounded public-price commitment, then solve location deviations and two-stage entry/public establishment. Re-do the prior-art absorption checks against Cremer, Marchand & Thisse (1991), Fernández-Ruiz (2020), Zhang & Li (2013), and mixed public-price studies; price continuation multiplicity is unlikely to be novel in itself.

Key background (abstracts / methodology):
- Cremer, Marchand, Thisse (1991): https://doi.org/10.1016/0167-7187(91)90004-5
- Fernández-Ruiz (2020): https://doi.org/10.1007/s12076-020-00249-y
- Zhang & Li (2013): https://pmc.ncbi.nlm.nih.gov/articles/PMC7080093/
- Kim (2022): https://doi.org/10.1080/1331677X.2021.1962385

Reproduction: python hotelling_entry_experimentation/code/stage5_price_subgame_audit.py
