# Stage 4 Verification Results

Date: 2026-09-01

Scripts:

- `research_reboot/stage4/verification/stage4_symbolic_verification.py`
- `research_reboot/stage4/verification/stage4_numerical_sanity.py`

Status: **PASS for the algebra of the minimal model.**

This is not a novelty verdict.

## 1. Symbolic identities verified

With

`D=1+4v(1-alpha)`

and

`n(alpha)=[sqrt(D)-1]/2`,

SymPy verifies exactly:

`n(1+n)=v(1-alpha)`.

Participation derivatives:

`dn/dalpha=-v/sqrt(D)`

and

`d2n/dalpha2=-2v^2/D^(3/2)`.

Regional objective:

`U=A alpha+B n(alpha)`.

The script verifies

`U'=A-Bv/sqrt(D)`

and

`U''=-2Bv^2/D^(3/2)`.

The candidate regional interior solution is

`alpha_R=1+1/(4v)-B^2 v/(4A^2)`

with implied participation

`n_R=(vB-A)/(2A)`.

Both the FOC and participation identity are verified exactly.

## 2. Planner identities verified

Outsider equilibrium net surplus is

`S_X=n^2/2`.

Planner welfare:

`W=U+n^2/2`.

The script verifies the decomposition

`W'=U'+n n'`

and the closed form

`W'=A-v/2-v(B-1/2)/sqrt(D)`.

It also verifies

`W''=v^2(1-2B)/D^(3/2)`.

For `B>1/2`, planner welfare is strictly concave.

The candidate planner solution is

`alpha_SP=1+1/(4v)-v(B-1/2)^2/[4(A-v/2)^2]`

with

`n_SP=(vB-A)/(2A-v)`.

The participation identity at this candidate is verified exactly. The planner FOC is verified using the positive-root relation valid on the target region `A>v/2`, `B>1/2`.

## 3. Joint interior region

A sufficient open parameter region is

`Omega={0<v<2, B>1/2, v/2+v(B-1/2)/sqrt(1+4v)<A<vB}`.

On this region:

- `0<alpha_SP<alpha_R<1`;
- `0<n_R<n_SP<1`;
- regional and planner optima are unique and interior;
- `U''<0` and `W''<0`.

The ordering follows especially cleanly from

`W'(alpha_R)=n(alpha_R)n'(alpha_R)<0`.

Since `W` is strictly concave, `alpha_SP<alpha_R`.

## 4. Local-welfare reversal

On `Omega`,

`U'(0)>0`,

`U'(1)<0`,

and `U''<0`.

Therefore there is a unique interior `alpha_R`, and for every `alpha>alpha_R`, increasing home bias reduces the regional sponsor's own payoff.

This is not a planner-weight artifact.

## 5. Fixed-participation mechanism test

With `n=nbar` fixed,

`U_F=A alpha+B nbar`,

so

`dU_F/dalpha=A>0`.

The regional sponsor chooses the maximal bias `alpha=1`; the local-welfare reversal disappears.

The planner may still dislike bias because outsider utility directly falls with `alpha`:

`W_F'=A-nbar*v/(1+nbar)`.

Hence a regional/planner wedge can survive without endogenous participation as a standard welfare-weight effect. This means the composition wedge by itself is not the distinctive result; only the sponsor-own-welfare reversal passes the mechanism-essentiality test.

## 6. Numerical sanity audit

Benchmark:

`v=1, B=1, A=0.8`.

Results:

- `n_R=0.125`
- `n_SP=0.333333...`
- `alpha_R=0.859375`
- `alpha_SP=0.555555...`
- `U_R(alpha_R)=0.8125`
- `U_R(1)=0.8`
- `W(alpha_SP)=0.833333...`
- `W(alpha_R)=0.8203125`

The numerical script also evaluates 25 admissible grid points spanning

- `v in {0.2,0.5,1.0,1.5,1.9}`
- `B in {0.6,0.8,1.0,1.5,2.0}`

with `A` chosen at the midpoint of the corresponding open admissible interval.

All 25 points satisfy:

`0<alpha_SP<alpha_R<1`,

`0<n_R<n_SP<1`,

`U_R(alpha_R)>U_R(1)`,

and

`W(alpha_SP)>W(alpha_R)`.

Numerical calculations are used only as a sanity check; the results above have analytic proofs.

## 7. Proposition-level verification status

- P1 home-bias composition wedge: **TRUE algebraically, but not sufficient for novelty and not fully mechanism-essential because a planner/regional wedge can remain under fixed participation through outsider welfare weights.**
- P2 interior participation tipping: **FALSE in the minimal model.** Participation declines smoothly and reaches zero only at `alpha=1`.
- P3 local-welfare reversal: **TRUE and mechanism-essential.** It disappears when participation is fixed.
- P4 neutrality commitment: **NOT PURSUED.** The sponsor already chooses its own interior optimum; adding a commitment game would exceed the Stage 4 freeze and is not warranted after the novelty/isomorphism concerns.

## 8. Verification conclusion

The minimal M10 mechanism does not collapse mathematically. Its strongest result is the participation-essential local-welfare reversal.

The remaining gate is economic novelty. The algebra alone does not show that this result is distinct from biased-intermediary/platform participation models.