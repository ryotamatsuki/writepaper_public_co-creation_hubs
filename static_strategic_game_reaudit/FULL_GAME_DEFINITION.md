# Full Static Game Definition

## Players and primitives

There are two regions `i=1,2`, regional governments `G_i`, public hubs `H_i`, a metropolitan private intermediary `H_T`, and a continuum or finite population of project types `theta`.

Each intermediary has access to a partner opportunity set `S_k`, `k in {1,2,T}`. For the public hubs the screening/evaluation investment `x_i` operates on `S_i` and changes the expected productive match value

\[
M_i(\theta,x_i;S_i).
\]

The metropolitan technology is fixed in the static baseline:

\[
M_T(\theta;S_T).
\]

This notation does not yet select a functional form. Stage 3 must choose one that retains economically meaningful opportunity-set geometry.

## Strategies

Regional government `G_i` chooses

\[
x_i\ge 0.
\]

The metropolitan private intermediary observes `(x_1,x_2)` and chooses

\[
p_T\ge 0.
\]

Private profit is conceptually

\[
\Pi_T=(p_T-c_T)D_T,
\]

with `c_T=0` available as a baseline normalization.

## Firm/project choice

A project type chooses

\[
h^*(\theta)\in\arg\max\{M_1(\theta,x_1;S_1),M_2(\theta,x_2;S_2),M_T(\theta;S_T)-p_T,0\}.
\]

The outside option is essential because it permits genuine extensive-margin match creation.

## Timing

1. `G_1,G_2` choose `(x_1,x_2)` simultaneously.
2. `H_T` observes public investments and chooses `p_T`.
3. Projects choose `H_1,H_2,H_T,0`.
4. Productive matches and surplus incidence realize.

The intended solution concept is subgame-perfect equilibrium, subject to the exact Stage 3 microfoundation.

## Regional and social welfare

At the conceptual level regional welfare is

\[
W_i=\text{resident-firm real productive surplus}
+\text{resident-partner incidence}
-C_i(x_i),
\]

with partner-incidence terms retained only if needed and microfounded.

The social planner counts total real productive surplus across regions and the metropolitan sector minus real resource costs. Access-price payments are transfers at the aggregate social level although their incidence can matter for regional objectives and strategic incentives.

## Opportunity-set geometry

A permissible primitive representation can have common and unique components, e.g.

\[
S_i=S^C\cup S_i^U.
\]

A descriptive overlap index such as

\[
\omega=\nu(S_1\cap S_2)/\nu(S_1\cup S_2)
\]

may be derived, but it cannot replace the primitive matching structure or be inserted as an arbitrary competition coefficient.