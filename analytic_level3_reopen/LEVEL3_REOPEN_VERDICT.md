# Level-3 Reopen Verdict

## Final classification

**NO-GO LEVEL 3 / GO LEVEL 2**

This verdict is materially different from L3-0.

## What L3-0 got wrong

The obstacle is **not mathematical non-characterizability**.

After polynomialization and denominator-sign bookkeeping, the regular-interior
reversal set is the projection of a semialgebraic set. Tarski–Seidenberg implies
that a primitive necessary-and-sufficient semialgebraic characterization exists.

Therefore "the condition may not exist" is rejected.

## What L3-1 established

1. Symmetry-first coordinates `(t,q_T,p)` remove `x` exactly.
2. The network parameter is invertible on the participation manifold:
   `delta=q_T(q_T-A_T)/(2(t q_T-a))`.
3. The generic symmetric private FOC reduces to a degree-7-in-`t`, quartic-in-`p`
   non-boundary core.
4. Canonical G public equilibrium elimination factors strongly; the target branch
   lies on a degree-31 algebraic component.
5. Canonical G `t` and `p` roots are exactly Sturm-isolated.
6. A component Gröbner basis makes `p` linear over the degree-31 `t` field.
7. B3 FOC, own derivative and cross derivative all reduce to polynomial sign
   blocks linear in `q_T`.
8. Raw symbolic size falls sharply after equilibrium-ideal reduction.

## Why Level 3 still fails

A Level-3 success requires the complete B3/G reversal decision, not only exact
G root geometry.

Two obligations remain open:

- matched-price coupling: the algebraic G price must select the distinct B3
  public equilibrium root with exact branch isolation;
- G sign query: the complete reduced G public own and cross derivatives must be
  sign-determined at the selected algebraic equilibrium.

Direct symbolic G cross construction reaches roughly 390k raw operations and
more than 425k after common-denominator assembly before final ideal reduction.
The direct resultant coupling the degree-31 G price component to the B3 core did
not complete within the 60-second resource cap.

No complete primitive projection, no finite primitive N&S condition set, and no
completed L3-R root/sign package was obtained.

## Strong Level-2 status

The L3-0 small-`beta` theorem remains valid and is now supported by a clearer
algebraic architecture. It is still the strongest proved theorem: on continuing
regular/SOC symmetric branches, sufficiently small positive `beta` implies
`BR_i^{B3'}>0>BR_i^{G'}`.

## Mathematical obstacle versus expository obstacle

- mathematical non-characterizability: **NO**;
- exact decidability in principle: **YES**;
- exact primitive projection computed: **NO**;
- current computational obstacle: **YES**;
- human-readable single threshold: **not found**.

## Exact next step

Do not repeat broad L3-1 exploration. A follow-up should be a narrowly scoped
**L3-2 component-wise exact sign engine**:

1. keep the degree-31 canonical/general component factored;
2. compute the G reduced Hessian modulo the augmented equilibrium ideal before
   common-denominator expansion, preferably with modular/subresultant arithmetic;
3. derive a linear/rational representation of `p_G` on each admissible component;
4. couple B3 using component-wise Gröbner/subresultant elimination rather than a
   monolithic resultant;
5. implement Thom/Sturm sign queries for `M_G`, `M_B3`, and both SOCs;
6. only then reopen primitive CAD or parametric `beta` boundary analysis.

A CAS with Singular/Sage/Mathematica/Maple/Redlog capabilities would materially
improve this next step.
