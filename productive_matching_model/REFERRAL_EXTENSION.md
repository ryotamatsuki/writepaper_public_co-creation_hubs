# Referral Extension M1

Let

`k=(1-p)^2 r`

be the incremental productive-match probability created by the foreign-specialist referral pool when both hubs are open. Since `h=2p-p^2`, feasibility gives `0<=k<=(1-p)^2` and `h+k<=1`.

The exact cross-entry effect is

`Gamma(k) = Gamma_0 + Delta k[8A+Delta(11-14h-7k)]/18`,

with

`Gamma_0 = -Delta p[8A+Delta(11-14p^2+7p^3)]/18 < 0`.

## Referral monotonicity

Differentiate with respect to `k`:

`d Gamma/dk = Delta[8A+Delta(11-14h-14k)]/18`.

Because `h+k<=1`,

`11-14h-14k >= -3`.

Hence under `A>Delta>0`,

`d Gamma/dk > Delta(8A-3Delta)/18 > 0`.

Thus referral strength strictly raises the incentive complementarity of public hub provision throughout the admissible domain.

## Unique threshold

If `Gamma((1-p)^2)>0`, there is a unique `k* in (0,(1-p)^2)` satisfying `Gamma(k*)=0`.

Let

`B = 8A+Delta(11-14h)`

and

`C = 8A+Delta(11-14p^2+7p^3)`.

Then `k*` is the smaller feasible root of

`7 Delta k^2 - B k + p C = 0`,

namely

`k* = [B-sqrt(B^2-28 Delta p C)]/(14 Delta)`.

The referral threshold is

`r*=k*/(1-p)^2`.

If the maximal-feasible referral value does not make `Gamma` positive, hubs remain substitutes for all `r`.

## Maximal-referral test

At `r=1`,

`Gamma_max = Delta[8A(1-3p+p^2)+Delta(4-33p+39p^2-14p^3)]/18`.

Thus a sign reversal exists exactly when this expression is positive.

## Exact interior witnesses

Take

`A=1, Delta=1/2, p=1/10`.

Then `h=19/100`.

### Negative open-set witness

At `r=1/10`, `k=81/1000`:

`D_0 = 1283/18000 ≈ 0.0712778`,

`D_1 = 485657/8000000 ≈ 0.0607071`,

`Gamma = -761087/72000000 ≈ -0.0105707`.

### Positive open-set witness

At `r=1/5`, `k=81/500`:

`D_0 = 1283/18000 ≈ 0.0712778`,

`D_1 = 2693/31250 ≈ 0.086176`,

`Gamma = 33521/2250000 ≈ 0.0148982`.

Both parameter points are strictly interior. Joint ±5% perturbations of all four primitives `(A,Delta,p,r)` preserve the negative sign in 100,000/100,000 valid draws around the first point and the positive sign in 100,000/100,000 valid draws around the second.

## Verdict

**Shared-pool/referral channel: ESSENTIAL for complementarity in this minimal model.**

With `r=0`, complementarity is impossible. With sufficiently strong structurally microfounded referral, complementarity occurs on nonempty interior open sets.