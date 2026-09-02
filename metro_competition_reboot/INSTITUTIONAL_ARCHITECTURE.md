# Institutional Architecture

## Players

- Regions `i in {1,2}` with regional governments `G_i`.
- Immobile local projects/firms with heterogeneous project need `theta`.
- Regionally anchored partner ecosystems `S_i`.
- Public matching intermediaries `H_i`, financed by `G_i`.
- One metropolitan/national private matching intermediary `H_T` with partner set/technology `S_T`.

## Core productive primitive

For a project `theta` and partner `s`, let `v(theta,s)` be the real productive surplus from collaboration. Values can be nonpositive. A hub does not add a separable amenity. It improves the process for identifying/evaluating feasible partners and selecting a productive collaboration from its accessible ecosystem.

A local ecosystem is a set or measure space `S_i`; its common component with the rival is primitive. A baseline overlap index may be reported as

`omega = nu(S_1 intersect S_2) / nu(S_1 union S_2)`,

or as a project-weighted analog. The model must be derived from the underlying sets/distributions; `omega` cannot appear as an arbitrary `-gamma x_i x_j` term.

## Local matching capacity

`x_i >= 0` is effective screening/evaluation capacity. It determines how many/effectively how deeply candidate partners from `S_i` can be evaluated for a project. If `Q_i(theta,x_i) subseteq S_i` denotes the evaluated opportunity set, the induced productive object is

`M_i(theta,x_i) = E[max{0, v(theta,s): s in Q_i(theta,x_i)}]`.

The expression is a contract for Stage 1, not yet a frozen stochastic specification.

## Metropolitan intermediary

Baseline `H_T` has exogenous matching technology `M_T(theta)` and chooses access price `p_T`. Its breadth must come from the capability distribution in `S_T`, not a free-standing 'Tokyo quality' constant.

## Firms do not move

Using `H_j` or `H_T` never relocates a beneficiary firm. This preserves a regional-incidence question distinct from firm-location competition.