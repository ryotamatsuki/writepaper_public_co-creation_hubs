# Stage 7 — Welfare Derivation

Project: `ryotamatsuki/writepaper_public_co-creation_hubs`

Workflow SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`

Stage 6 input SHA: `34a955faa1814ac8312051b971c1f75d59aad1ee`

## 1. Canonical notation

Let

\[
d\equiv A_L-A_M,\qquad q\equiv 1-2\beta.
\]

The canonical interior allocation is

\[
h=\frac{d-2\beta-\beta\tau}{q},\qquad x=h-\tau,\qquad s=\frac{d-2\beta}{q}.
\]

The interior domain is

\[
2\beta+\tau(1-\beta)<d<1.
\]

Core payoffs are

\[
C_0=A_M+2\beta,
\]

\[
C_1=A_M+\beta(2-2h+\tau),
\]

\[
C_2=A_M+\beta(2-2s).
\]

There is no separate consumer population in the canonical model. Hence **consumer surplus as a separate welfare object is not applicable**. Aggregate welfare is the realized utility of the two unit masses of firms/users, net of real establishment costs.

## 2. Aggregate welfare from individual utilities

### No peripheral entry

All firms use the core, so

\[
B_0\equiv SW_{00}(F=0)=2C_0.
\]

Thus

\[
SW_{00}=2C_0.
\]

### One peripheral entry

With region 1's platform present, region-1 firms with \(c<h\) use their own platform, while region-2 firms with \(c<x=h-\tau\) use region 1's platform. Integrating the gain over the common core payoff gives

\[
B_1\equiv SW_{10}(F=0)=2C_1+\frac{h^2+x^2}{2}.
\]

Hence

\[
\boxed{SW_{10}=SW_{01}=2C_1+\frac{h^2+(h-\tau)^2}{2}-F.}
\]

### Two peripheral entries

Each region's firms with \(c<s\) use the local platform. Thus

\[
B_2\equiv SW_{11}(F=0)=2C_2+s^2,
\]

and

\[
\boxed{SW_{11}=2C_2+s^2-2F.}
\]

The previously reported Stage-5 welfare formulas are therefore verified from primitives.

## 3. Regional/private entry benefits

For a region considering its own establishment:

\[
D_0=\left(C_1+\frac{h^2}{2}\right)-C_0
=\beta(\tau-2h)+\frac{h^2}{2},
\]

when the other region has not entered, and

\[
D_1=\left(C_2+\frac{s^2}{2}\right)-\left(C_1+\frac{x^2}{2}\right)
=\beta(2h-2s-\tau)+\frac{s^2-x^2}{2},
\]

when the other region has entered.

The strategic interaction is

\[
\Gamma=D_1-D_0.
\]

## 4. Cross-region welfare effects

Define the effect of region 1's entry on region 2 when region 2 has not entered:

\[
X_0\equiv \left(C_1+\frac{x^2}{2}\right)-C_0
=\beta(\tau-2h)+\frac{x^2}{2}.
\]

Define the effect of the second entry on the region that already has a platform:

\[
X_1\equiv \left(C_2+\frac{s^2}{2}\right)-\left(C_1+\frac{h^2}{2}\right).
\]

Symbolic simplification yields

\[
\boxed{
X_1=-\frac{\beta\tau\,[2-2d+\beta\tau]}{2(1-2\beta)^2}<0
}
\]

throughout the canonical interior domain because \(\beta>0\), \(\tau>0\), and \(d<1\).

This corrects the sign typo in the Stage-7 prompt: the candidate expression requires the leading minus sign.

The first-entry effect \(X_0\) is sign-ambiguous. A first peripheral platform may benefit firms in the other region through cross-use, while also reducing the network value of the common core.

## 5. Social marginal establishment benefits

Define

\[
G_0\equiv B_1-B_0,
\qquad
G_1\equiv B_2-B_1.
\]

Then exactly

\[
\boxed{G_0=D_0+X_0,}
\]

\[
\boxed{G_1=D_1+X_1.}
\]

Hence the private/social wedges are precisely the effects of one jurisdiction's establishment on firms originating in the other jurisdiction.

In cutoff form,

\[
G_0=2\beta(\tau-2h)+\frac{h^2+x^2}{2},
\]

and

\[
G_1=2\beta(2h-2s-\tau)+s^2-\frac{h^2+x^2}{2}.
\]

## 6. Key identities

A direct rearrangement of definitions gives

\[
\boxed{X_1-X_0=\Gamma.}
\]

Therefore

\[
\boxed{G_1-G_0=2\Gamma.}
\]

This is the central Stage-7 welfare identity: the same threshold \(A^*\) that reverses private bilateral entry interaction also reverses the ordering of social marginal entry benefits, with twice the signed difference.

Let

\[
\bar G\equiv\frac{G_0+G_1}{2}.
\]

Using the identities above,

\[
\boxed{\bar G=D_0+X_1<D_0.}
\]

This strict inequality is crucial for the welfare characterization of the coordination region.

## 7. Planner problem

Net welfare for \(N=0,1,2\) peripheral platforms is

\[
V_0=B_0,
\qquad
V_1=B_1-F,
\qquad
V_2=B_2-2F.
\]

### Strategic-substitutes region: \(\Gamma<0\)

Because \(G_1<G_0\), social marginal benefits are decreasing. Subject to \(F>0\):

- \(N^{SP}=2\) if \(F<G_1\) (when this interval is positive);
- \(N^{SP}=1\) if \(G_1<F<G_0\);
- \(N^{SP}=0\) if \(F>G_0\).

### Strategic-complements region: \(\Gamma>0\)

Because \(G_1>G_0\), one-entry is never globally optimal. Comparing the endpoints gives

\[
V_2-V_0=G_0+G_1-2F=2(\bar G-F).
\]

Hence

- \(N^{SP}=2\) if \(F<\bar G\) (when \(\bar G>0\));
- \(N^{SP}=0\) if \(F>\bar G\);
- \(N^{SP}=1\) is never the strict social optimum.

Because

\[
\bar G<D_0<D_1
\]

in the coordination region, the social full-entry threshold lies strictly below the private threshold that makes no entry a best response.

## 8. Welfare propositions

### W1 — Efficient anti-coordination exists

The private one-entry equilibrium region is

\[
\max\{0,D_1\}<F<D_0.
\]

The planner chooses one entry when

\[
\max\{0,G_1\}<F<G_0.
\]

Their overlap is nonempty on an open parameter set. For example, with

\[
\beta=0.05,\quad \tau=0.20,\quad A_L-A_M=0.60,
\]

we obtain

\[
D_0=0.1037654,\quad D_1=0.0838889,
\]

\[
G_0=0.1186420,\quad G_1=0.0788889.
\]

Thus \(F=0.09\) yields the asymmetric one-entry Nash equilibria and \(N^{SP}=1\). All inequalities are strict, so this is an open-set existence result.

### W2 — Coordination can generate socially excessive proliferation

Suppose \(\Gamma>0\) and \(D_1>0\). Whenever

\[
\max\{0,D_0\}<F<D_1,
\]

both \((0,0)\) and \((1,1)\) are Nash equilibria. But \(F>D_0>\bar G\), so

\[
SW_{00}>SW_{11}.
\]

Hence the full-entry coordination equilibrium is socially excessive. This is a theorem, not a numerical example.

A stronger result holds for

\[
\max\{0,\bar G\}<F<D_0
\]

when \(D_0>0\): full entry is the unique Nash equilibrium but the planner chooses no peripheral entry.

### W3 — No full under-entry coordination failure

There is no parameter region in which \((0,0)\) is a Nash equilibrium while \(N^{SP}=2\).

- If \(\Gamma<0\), no entry requires \(F>D_0>D_1>G_1\), which excludes social full entry.
- If \(\Gamma>0\), no entry requires \(F>D_0>\bar G\), which again excludes social full entry.

Thus the proposed W3 result is **FAIL / impossible in the canonical model**.

### Additional welfare corollary — under-provision of the first platform can occur

In the substitutes region, if \(X_0>0\), then \(G_0>D_0\). For

\[
D_0<F<G_0,
\]

no entry is the Nash equilibrium while the planner chooses one peripheral platform. Thus under-entry can occur only at the first-entry margin; the canonical model does not produce a no-entry equilibrium when full entry is socially optimal.

## 9. Social game-form transition

Because

\[
\boxed{G_1-G_0=2\Gamma,}
\]

the threshold \(A^*\) is not merely a private strategic threshold. At exactly the same \(A^*\), social marginal establishment benefits change from decreasing to increasing:

\[
A_M<A^*\Rightarrow G_1<G_0,
\]

\[
A_M>A^*\Rightarrow G_1>G_0.
\]

The welfare and strategic thresholds therefore coincide in sign, although the private and social entry-cost cutoffs do not coincide.

## 10. Verification

SymPy independently verified:

- all cutoff substitutions;
- \(SW_{00},SW_{10},SW_{11}\);
- \(G_0-D_0=X_0\);
- \(G_1-D_1=X_1\);
- the closed form and strict sign of \(X_1\);
- \(X_1-X_0=\Gamma\);
- \(G_1-G_0=2\Gamma\);
- \((G_0+G_1)/2=D_0+X_1\).

A 200,000-draw admissible-parameter stress test found:

- \(X_1<0\) in 200,000/200,000 draws;
- numerical error in \(G_1-G_0=2\Gamma\) below \(1.2\times10^{-16}\);
- zero cases in which no entry could be a Nash equilibrium while full entry was socially optimal;
- every sampled coordination-multiplicity case had \(N^{SP}=0\).

Monte Carlo frequencies are verification only, not economic probability statements.
