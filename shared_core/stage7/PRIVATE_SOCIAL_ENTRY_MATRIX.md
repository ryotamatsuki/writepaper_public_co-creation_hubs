# Stage 7 — Private vs Social Entry Matrix

## Definitions

Private/regional gross entry incentives:

\[
D_0=W_i(1,0;F=0)-W_i(0,0;F=0),
\]

\[
D_1=W_i(1,1;F=0)-W_i(0,1;F=0).
\]

Cross-region effects:

\[
X_0=W_j(0,1;F=0)-W_j(0,0;F=0),
\]

\[
X_1=W_j(1,1;F=0)-W_j(1,0;F=0).
\]

Social gross marginal benefits:

\[
G_0=D_0+X_0,
\qquad
G_1=D_1+X_1.
\]

The canonical identities are

\[
X_1<0,
\qquad
X_1-X_0=\Gamma,
\qquad
G_1-G_0=2\Gamma.
\]

## Strategic-substitutes region: \(\Gamma<0\)

Private thresholds satisfy

\[
D_1<D_0.
\]

Social thresholds satisfy

\[
G_1<G_0.
\]

Because \(X_1<0\),

\[
G_1<D_1.
\]

For positive fixed cost \(F\):

| Fixed-cost region | Nash outcome | Social optimum | Interpretation |
|---|---|---|---|
| \(0<F<\min\{G_1,D_1\}\), if positive | \(N^{NE}=2\) | \(N^{SP}=2\) | efficient full entry |
| \(G_1<F<D_1\), where nonempty | \(N^{NE}=2\) | \(N^{SP}\le1\) | excessive second entry |
| \(\max\{0,D_1\}<F<D_0\) | asymmetric one-entry NE | \(N^{SP}=1\) if \(F<G_0\), otherwise 0 | efficient specialization or excessive first entry |
| \(D_0<F<G_0\), possible iff \(X_0>0\) | \(N^{NE}=0\) | \(N^{SP}=1\) | under-provision of first entry |
| \(F>\max\{D_0,G_0\}\) | \(N^{NE}=0\) | \(N^{SP}=0\) | efficient no entry |

The anti-coordination region therefore permits both efficient one-platform specialization and, depending on the sign of the first-entry cross-region effect \(X_0\), under- or over-provision at the first-entry margin.

## Strategic-complements region: \(\Gamma>0\)

Private thresholds satisfy

\[
D_0<D_1.
\]

Social marginal benefits satisfy

\[
G_0<G_1.
\]

Define

\[
\bar G=\frac{G_0+G_1}{2}.
\]

The canonical welfare identity implies

\[
\boxed{\bar G=D_0+X_1<D_0<D_1.}
\]

The planner never strictly chooses exactly one platform in this region. For positive \(F\):

| Fixed-cost region | Nash outcome | Social optimum | Interpretation |
|---|---|---|---|
| \(0<F<\bar G\), if \(\bar G>0\) | full entry unique | full entry | efficient full entry |
| \(\max\{0,\bar G\}<F<D_0\), if \(D_0>0\) | full entry unique | no entry | severe excessive proliferation |
| \(\max\{0,D_0\}<F<D_1\) | \((0,0)\) and \((1,1)\) | no entry | no-entry equilibrium efficient; full-entry equilibrium excessive |
| \(F>D_1\) | no entry unique | no entry | efficient no entry |

If \(D_0\le0<D_1\), the coordination-multiplicity region begins at \(F=0\); because \(\bar G<D_0\le0\), the social optimum is no entry for every positive \(F\) in that multiplicity region.

## Important negative result

There is **no** open region in the canonical model in which no entry is a Nash equilibrium while full entry is socially optimal.

Thus the coordination transition does not create the familiar under-entry coordination failure. Instead, the network-dependent residual core produces a one-sided welfare problem: self-reinforcing peripheral entry can generate excessive proliferation.

## Economic interpretation

Strategic complementarity does not mean that one region's entry creates a positive welfare spillover on the other. In the coordination region both cross-region effects are negative:

\[
X_1<0,
\]

and since \(X_0=X_1-\Gamma\),

\[
\Gamma>0\Rightarrow X_0<0.
\]

The complementarity is in **entry incentives conditional on the other's entry state**, not in direct cross-region welfare. This distinction is important for exposition and prevents an incorrect interpretation of the result as cooperative spillovers.
