# Stage 11R P1j — model-derived economic insight: technology makes second-best public access more valuable despite adverse local innovation

**Research date:** 2026-10-28 JST. Baseline is JICT Draft PR #61 (head `cb92ec3cc5bb6fb94b20420bc5085605064d00cc`). Sources are the unchanged `04_propositions.tex`, `04_diagnostic.tex` and existing Stage33 Lean partial formalization. This independent stage **does not alter** the model, main manuscript, parameter calibrations or Lean. Parent chain PR #62–#66 already separates pricing/game-theoretic novelty from model-specific contributions.

## 1. Executive verdict

A more theoretically robust and economically intelligible result than the narrow three-scenario entry switch was available **inside the already frozen model**:

> **Technology-induced screening/access inversion:** Even when `K_i=M_i+omega_i*Pi<0` (successfully commercializing a venture reduces region i's net venture-related surplus), a rise in skill complementarity `lambda` can **strictly increase** the maximum fixed cost `T_i(lambda)` that region i would rationally pay to add free basic service alongside a fee-setting private specialist.

The channel is **endogenous private fee and founder exclusion**, not positive regional commercialization. The government restores some excluded basic-service value, while the additional matched ventures can be locally adverse. Both the baseline technology term and the skill-success increment are negative over an explicitly verified **broad existing-primitive continuous rectangle**, yet the policy threshold strictly rises.

This does NOT say government first-stage entry becomes strictly optimal for every fixed government cost (the final binary switch at baseline costs is narrow); it is a **conditional policy-value comparative static**, with the other region full (G/B) and the local private supplier P versus B. This finding already appears qualitatively in PR #61's `diagnostic-derivative`; **the incremental contribution of P1j is a clean primitive inequality, its nonempty economically interpretable incidence band, and a strong parameter-continuum rational certificate**, not discovery of a new generic platform or public-good principle.

## 2. Exact primitive theorem under weak bilateral matching

Maintain independent uniform founder skill/valuation `theta_i`, `0<u_i<d_i=Delta_i`, positive downstream venture-profit `Pi`, nonnegative `R0`, positive `N_i,N_j`, partner `j` with full participation `(m_j,q_j)=(1,1/2)`, nonnegative `gamma=rho*mu` and a fixed `rho` as `gamma->0+`. Keep `lambda` in any fixed compact interval in the interior success-probability domain `p_b+lambda<=1`.

Put:
```
C=gamma*N_i*N_j,
U_i=u_i+gamma*N_j*omega_i*R0,
eta_i=gamma*N_j*omega_i*Pi/2,
D_i=Delta_i+eta_i*lambda,
ell_i=(D_i-U_i)/(2D_i),
K_i=M_i+omega_i*Pi,
S_i0=p_b*M_i+omega_i*R0.
```

The **frozen-model** public addition threshold from P to B is
```
T_i(lambda)=N_i[u_i*ell_i-Delta_i/8+Delta_i*ell_i^2/2]
             +C[S_i0*ell_i+lambda*K_i*ell_i^2/4].
```

Differentiate and use `ell'_i=U_i*eta_i/[2*D_i^2]`. Exact:
```
T'_i(lambda)=
   [N_i(u_i+Delta_i*ell_i)+C*S_i0]*ell'_i
   +C*K_i*(ell_i^2+2*lambda*ell_i*ell'_i)/4.
```

For `gamma->0+` at **fixed rho** and any fixed compact admissible lambda range (with sufficient positive distance from participation and success corners), uniformly:
```
T'_i(lambda) =
 gamma*N_i*N_j/(16*Delta_i^2)
 * [ M_i*(Delta_i-u_i)^2 + omega_i*Pi*(Delta_i^2+3*u_i^2) ]
 + O(gamma^2).
```

Sympy algebra was checked independently: substituting `D=Delta_i+eta*lambda` in the exact polynomial identity has ZERO symbolic remainder, and the small-`gamma` expression reduces exactly to the bracket above.

Hence, if the bracket is strictly positive (negative), for all sufficiently small positive `gamma`, the public access threshold increases (decreases) with lambda uniformly over the compact lambda domain, not just at a selected lambda.

When `M_i<0`, the following primitive **open profit-incidence band** yields simultaneously `K_i<0` and `T'_i(lambda)>0`:
```
(-M_i)*(Delta_i-u_i)^2/[Pi*(Delta_i^2+3*u_i^2)]
    < omega_i < (-M_i)/Pi.
```
Also intersect with `0<=omega_i<=1`. The algebraic interval is nonempty if `u_i>0` because
`(Delta_i-u_i)^2 < Delta_i^2+3*u_i^2`. A *feasible share* interval requires the lower bound to be <1. This is an economically meaningful new **model-specific primitive condition**, not a general normal-form game theorem.

Define `r_i=u_i/Delta_i`. The lower critical appropriation share is
```
omega_i^crit = (-M_i/Pi)*(1-r_i)^2/(1+3*r_i^2).
```
For `0<r_i<1`, its derivative w.r.t. `r_i` is strictly negative:
```
d omega_i^crit / dr_i = (-M_i/Pi)*[-2*(1-r_i)*(1+3*r_i)]/(1+3*r_i^2)^2 < 0.
```
**Interpretation:** where access to the basic service itself matters more, less local venture-profit appropriation is needed for technical innovation to intensify private screening sufficiently to make public access more attractive. At `omega_i=0`, the price-selection channel disappears, and with `K_i=M_i<0`, the government value of lambda **falls** on any interior branch (an exact statement, not only a weak-gamma limit).

## 3. Stronger exact finite-gamma parameter-continuum theorem

Using original PR #61 exact Cournot numbers `M_A≈-0.01579743681025276`, `M_B≈-0.03077500469131169`, `Pi≈0.11066084465254933`, `p0=0.05`, `u_B=0.28`, `gamma=0.1`, and partner A fully covered, the weak-gamma critical share interval is
```
omega_B^crit ≈ 0.11671639279397
omega_B^upper = -M_B/Pi ≈ 0.27810202233625
```

**New exact rational inclusion proof (`incidence_screening_sign_certificate.py`):** for **every single point**, not just sampled points, in
```
0.13 <= omega_B <= 0.27,  omega_A=1-omega_B,
0 <= lambda <= 0.5,
gamma=0.1, original other economics,
A founder participation (m_A,q_A)=(1,1/2),
```
all of the following are simultaneously true:
1. The local private-only price optimum is strictly interior `0<ell_B<1/2`.
2. Region B's market-plus-incumbent baseline innovation incidence `S_B0=p_b*M_B+omega_B*R0` is **strictly negative**.
3. Region B's extra successful-project skill impact `K_B=M_B+omega_B*Pi` is **strictly negative**.
4. `dT_B/dlambda>0` **strictly throughout the entire two-dimensional rectangle**.
5. The national full-coverage innovation income incidence `M_A+M_B+Pi>0` is independent of omega_B and positive at the original Cournot values.

The rational checker partitions omega into 10 closed subintervals and lambda into 10 closed subintervals. Each cell's exact Fraction interval extension includes all points in that subrectangle, and asserts that all 100 lower-bound derivative numerators are strictly positive; local K and S0 remain strictly negative, and all interior cutoff bounds strict. Unlike a Monte Carlo or dense point grid, this covers the entire closed parameter rectangle. Positive derivative lower bound `>0.0000334036961` in the monetary model units for **all** cells, based on exact outward inclusion; maximum (least negative) K_B `<-0.0008965766`; max S0 `<-0.0001543539`. These are conservative certified bounds.

The exact derivative numerator is the elementary polynomial identity
```
(16*D_i^3/C)*T'_i =
 omega_i*Pi*U_i * [4*D_i*u_i + 2*Delta_i*(D_i-U_i)]
 +8*D_i*U_i*eta_i*S_i0
 +K_i*(D_i-U_i) *
     [D_i*(D_i+U_i)-2*Delta_i*U_i],
```
with `D_i=Delta_i+eta_i*lambda`; equality follows by differentiating the frozen closed-form threshold. All coefficients and interval arithmetic are exact Fraction values. This is an **analytic derivative identity + proof-oriented exact interval certificate**, not merely a comparison at 100 sample points. The checker also verifies identity as exact equality at the three original-lambda rational positions.

At original `omega_B=.2`:
```
lambda=0:   T_B≈0.040425071007120596, T'_B≈0.00021374591618160668
lambda=0.25:T_B≈0.04047843897148336,  T'_B≈0.00021319799482516329
lambda=0.5: T_B≈0.04053167010244885,  T'_B≈0.00021265124845710635
```
The *discrete baseline public entry switch* at unchanged F_GB=.0405 is still only supported by a narrow cost range because its threshold moves by just ~0.0001066 over lambda 0->0.5. **Broad primitive comparative-static validity must never be confused with broad discrete entry-switch robustness.**

## 4. Additional deductions for interpretation and tests

### (a) A three-sign incidence paradox

At original `omega_B=.2`, the local fully covered marginal technical success welfare sign is `K_B=M_B+omega_B Pi≈-0.008642836<0`, whereas the national successful-product incidence is `M_A+M_B+Pi≈+0.0640884>0`. At the same time `dT_B/dlambda>0` for lambda in [0,0.5] by the above exact continuous certificate.

This simultaneously establishes:
1. **More valuable technical skill complementarity benefits the country** at given full participation,
2. **Its direct commercialisation effect harms region B** at given full participation, but
3. **The same technical advance increases B government's willingness to provide access**, because the indirect private screening component is stronger.

The three signs are not logically contradictory: #1–2 concern project profits/displacement in the *product market*, while #3 compares the *public and privately rationed participation regime with endogenous posted fees*. Proof of each belongs to a different conditional object. No universal qualitative prediction is implied.

### (b) The marginal restored pair is lower skill, not worthless

Under a given private-only cutoff `ell_i>0` and full partner participation, free entry restores the low-type group `theta_i in [0,ell_i)` with **mean skill `ell_i/2`**. The partner's mean is `1/2`, so the additional match's **mean skill product** is `ell_i/4`, which is lower than the full-population mean `1/4`. This is the microfoundation of
```
additional match mass = C*ell_i;
extra ability-weighted matches = C*ell_i^2/4;
extra expected successful projects = C*ell_i*(p_b+lambda*ell_i/4) > 0.
```
Therefore under the conditional P->B intervention **innovation successes increase even while B's regional project-market welfare declines** if S_B0<0,K_B<0. The project success rate for the *marginal* added pairs is lower than the overall full-participation rate; this is not an economic judgment that individual lower-skilled firms are unproductive.

### (c) Private fees are transfers; access exclusion—not retained rents—is the welfare channel

The founder's fee payment is an intraregional transfer to the intermediary and disappears from regional `W_i`. Thus the policy-value increase **cannot** be interpreted as a pure transfer from private rents to taxpayers. It instead corrects distorted *real service participation* that arises from monopoly fee maximisation, though co-provision may shift some customers away from the specialist premium and thereby lower real specialist-service surplus. The model only compares free basic public service to no public basic service, not to optimal rate regulation, vouchers, subsidies or investment in a physical hub.

### (d) Why an unobserved public opening is not always predicted

Under fixed original public cost `F_GB=0.0405`, a public entry switch occurs only after a threshold `lambda≈0.351184` numerically; that crossing is not certified over every `lambda` for full four-player equilibrium (the conditional `T_B(lambda)` rises continuously, but full off-path P/P fees for all lambda have not been uniformly certified). Thus the robust result is **policy-value monotonicity**, not an all-lambda **equilibrium-entry** theorem. The original baseline three cases are independently certified strict SPNE. Threshold widths are economically narrow and the Beta(2,2) ability robustness check is negative.

## 5. Adjacent original literatures: exact economic objects and nonabsorption limits

| Closest original paper | Players/choice and its actual screening/technology mechanism | What our result does beyond it; what is already known |
|---|---|---|
| [Hagiu (2011), *Quantity vs. Quality: Exclusion by Platforms with Network Effects*](https://www.hbs.edu/ris/publication%20files/11-125.pdf) | A profit-max platform selects **price AND an explicit exclusion/quality admission threshold**. Users care about network size and **average participant quality**, and the optimal quality threshold may rise at the cost of quantity. | **Selection and quality–quantity welfare conflict are known.** Our provider **does not choose a separate quality exclusion rule**; fee endogenously excludes types, and new R&D pair products enter two downstream Cournot markets with locally different signs. Hagiu does not provide the exact threshold `M(Delta-u)^2+omegaPi(Delta²+3u²)`. This absence is not a general nonembedding proof. |
| [Peitz & Sato (2026), *Asymmetric Platform Oligopoly*](https://doi.org/10.1111/1756-2171.70008) | Competing two-sided **commercial** platforms set membership prices for user groups facing logit taste and logarithmic within/cross-group network returns, with possible **free outside option and partial participation**. | General fee-participation feedback, exclusion, entry and welfare differences **already known**. Our intervention is independent welfare-maximising government's zero-priced basic service and region-specific commercialisation/Cournot losses. Their outside option can itself be freely supplied but platform competition vs co-provision is not mathematically equivalent by label alone. |
| [Mekonnen (2026), *Luck Out or Outpay? Competing with a Public Option*](https://arxiv.org/abs/2507.12779) | A capacity-constrained free public provider competes with a private monopolist who withholds supply; congestion and rationing change willingness to pay and private price. | Public expansion leading to private exclusion / higher fees / welfare losses **already known**; our mechanism lacks congestion/rationing and premium quality affecting match success. Instead, innovation skill complementarity strengthens private screening and the basic-service government's access option despite region's negative downstream venture incidence. |
| [Kalamov & Staal (2023), regional public goods](https://doi.org/10.1016/j.regsciurbeco.2023.103917) | Local governments' public-good choices interact via cross-region externality curvature and fiscal consequences. | Cross-jurisdiction positive/negative welfare wedges and strategic substitutes/complements are known. Our `T'_i>0` despite `K_i<0` is due to **endogenous private fee exclusion**, a distinct *model-specific implementation*, not a generic new principle. |
| [Zhang & Li (2013), *Endogenous R&D spillover and location choice in a mixed oligopoly*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7080093/) | Public/private **product-market firms** choose R&D, plant locations and prices. Public enterprises compete in products directly; R&D spillovers depend on location. | Mixed public/private R&D and region-specific spillovers already studied, but our government **does not compete as a public producer** and private premium service **does not raise R&D productivity directly**. The derived premium-price screening correction is different. |
| [Rossi et al. (2026), EDIH business models](https://doi.org/10.1093/icc/dtag042) | Empirical study of territorial/technical focus and mixed public/private innovation intermediary funding, not Cournot or Bayesian type-sorting theorem. | Institutional relevance for co-provision, not a theoretical priority certificate. |

**Editorial novelty verdict:** This **primitive-level, open-parameter sign condition** is a meaningful analytical strengthening of PR #61's already stated derivative formula and potentially deserves a *short corollary immediately below that formula* if main-paper length permits, but cannot be advertised as an entirely new platform/game-theoretic law. A referee might consider it a focused applied-theory contribution when tightly linked to the already verified original three-scenario outcomes; formal originality over arbitrary general spillover- and nonlinear-participation frameworks is **not certified**. Particularly, `q=phi(m)` so a general nonlinear network model may encode the same result.

**Manuscript policy:** Keep frozen PR #61 untouched in this exploratory PR. Do not add another headline proposition automatically. A useful revised central narrative would stress the *robust screening-induced comparative static*, not rely on the exceptionally narrow discrete entry-switch as if economically widespread. Remind readers that the reduced-form condition only holds in the maintained uniform skill==premium taste setup.

## 6. Reproducibility

Run in repository root:
```bash
python hotelling_entry_experimentation/production/government_game_exploration/incidence_screening_sign_certificate.py
```
It writes `incidence_screening_sign_results.json` (not committed, retained as a CI artifact). GitHub Actions workflow `.github/workflows/stage11r-p1j-incidence-screening.yml` runs exact interval proof + PR61 baseline global P/P pricing/entry + P1h finite-gamma all-four-player tests. No Lean verification of the NEW interval theorem, and no claim of external human referee signoff.
