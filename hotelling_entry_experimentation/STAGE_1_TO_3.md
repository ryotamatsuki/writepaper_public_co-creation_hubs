# Research log — Stage 1 to Stage 3 (as of 2026-10-27)

All claims below are archived from exploratory discussion. "Derivation" means a formula deduced under its specified restricted game, NOT a globally certified subgame-perfect equilibrium of the full model.

## Starting point
Quadratic Hotelling service-characteristic line x in [0,1], mismatch cost t(x-l_j)^2, public hub maximizes social welfare W, private hub maximizes profit. Location/characteristic choice then price competition. Goal: derive whether and how a public facility should be established, while allowing private entry.

### Stage 1 — Matching scope and investment
- Experimentation partners n_j in {1,2}; proposed private investment coefficient r_n = alpha*n/(1+sigma*(n-1)).
- Simplified optimizing choice: e_j = r_(n_j)*d_j with hub share d_j; private net investment contribution r_n^2*d_j^2/2.
- Candidate user service benefit: u_n(d)=eta*n*d + r_n^2*d^2/2 - k*n.
- Candidate social value: s_n(d)=(eta+chi)*n*d + B_n*d^2 - k*n, B_n=alpha*n*r_n-r_n^2/2.
- Earlier discussion reported an interior candidate in which the public hub chooses n_G=1, private hub n_P=2 (a about 0.2235, b about 0.6967), and nonmonotone scope selection as matching cost k varies.
- IMPORTANT: No executable Stage 1 optimizer, full deviation checks, price-equilibrium branch checks, or provenance of those numbers was provided to this repository. Treat every such number as UNVERIFIED. Do not cite Stage 1 numerical claims as a theorem or replicable equilibrium.

### Stage 2 — Government establishment and private entry
- Timing introduced: government build/no-build, private entry at fixed cost F_P, characteristics, prices, firm participation and investment.
- Fixed-location benchmark a=0, b=1, private service value v+delta, public value v. Full coverage and interior duopoly require explicit restrictions, with t>0, 0<=delta<t.
- With public welfare maximization and private profit maximization, price candidate p_G=p_P=t+delta, private market share (t+delta)/(2*t), pre-entry private profit M*(t+delta)^2/(2*t).
- Private-only monopoly under partial coverage v<3*t: p=2*v/3, private pre-entry profit (2*M*v/3)*sqrt(v/(3*t)).
- Private entry may be induced or deterred by the public service depending on private profitability in the relevant counterfactual.
- Local authority chooses build iff welfare under the build-induced private-entry equilibrium minus public fixed cost F_G is at least no-build welfare. The entry condition and social desirability of entry are separate objects.
- The generic entry-welfare discrepancy is not novelty (Mankiw–Whinston free-entry benchmark).

### Stage 3 — Endogenize spillover via pilot investment, *conditional on simplified timing and fixed positions*
- Local government chooses public establishment g in {0,1} and a pilot scope n in {1,2}. Participating pilot firm invests e, then a private hub decides entry, followed by hub price competition and firm service adoption.
- Pilot firm objective r_n*e-e^2/2; best reply e_n=r_n.
- Private hub spillover per participant: delta_n=gamma*n*e_n = gamma*alpha*n^2/(1+sigma*(n-1)).
- Direct pilot social surplus: H_n=r_n^2/2+B*n*r_n-k*n, EXCLUDING the downstream private spillover value to avoid double counting.
- In the restricted interior full-coverage public/private pricing game, private pre-entry profit Pi_P(n,M)=M*(t+delta_n)^2/(2*t).
- Private entry requires F_P <= Pi_P. Gross social gain from private entry compared with public-only in this benchmark is DeltaS=M*(t+delta_n)^2/(4*t).
- Under n=1 deterring entry and n=2 inducing it, with delta_2>delta_1 and H_2>H_1, a sufficient inequality for the public authority to prefer scope 1 is F_P>H_2-H_1+DeltaS(n=2,M). All comparisons are **within the simplified equilibrium regimes**.
- A candidate nonmonotone scope schedule 2->1->2 arises as market mass M crosses private entry thresholds M_n=2*t*F_P/(t+delta_n)^2, with additional welfare restrictions. This is a conditional fixed-location statement, NOT a fully endogenous-location theorem.

### Exact sample parameterization (see code)
t=1, v=2, alpha=0.30, sigma=1.40, gamma=1, B=0.50, k=0.04, F_P=1, F_G=0.60.
For n=1: e=0.300, delta=0.300, H=0.155, private pre-entry profit at M=1 is 0.845 -> NO private entry.
For n=2: e=0.250, delta=0.500, H=0.20125, private pre-entry profit at M=1 is 1.125 -> private ENTRY.
Net social welfare at M=1 (as computed in code):
- no public, private-only: ~0.451549
- public with n=1, public-only: ~1.221667
- public with n=2, public+private: ~0.830417
Private-entry mass thresholds: M_2 = 0.88888889 and M_1 = 1.18343195.
All amounts are normalized hypothetical model values, not empirical evidence.

## Claims NOT established
1. A fully endogenous location (a,b), scope (n_G,n_P), price, investment, entry and establishment subgame-perfect Nash equilibrium.
2. Global deviations and boundary equilibria in any joint location/price model.
3. Genuine endogenous bargaining/hold-up: r_n is a reduced-form investment recovery coefficient.
4. Whether physical hub construction dominates subsidy, contract, or digital alternatives.
5. Novelty against mixed oligopoly, spatial entry, two-sided platforms, relation-specific investment and over-entry literatures.

## Working research judgment
Use the fixed-location calculations as a unit-tested diagnostic benchmark. Never turn the illustrative over-entry result or 2->1->2 pattern into an unqualified publication claim.
