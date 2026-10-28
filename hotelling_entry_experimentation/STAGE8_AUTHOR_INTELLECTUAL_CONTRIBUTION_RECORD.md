# Stage 33 — Author intellectual-contribution and final theory acceptance record

Date prepared: 2026-10-28 JST
Record type: `AUTHOR_INTELLECTUAL_CONTRIBUTION_RECORD`
Status: `PARTIALLY DOCUMENTED — AUTHOR MATHEMATICAL REVIEW PENDING`
Theory scope ID: `STAGE33-2026-10-28-CANONICAL-SCOPE-RC1`
Canonical specification: `hotelling_entry_experimentation/STAGE8_CANONICAL_THEORY_SCOPE_LOCK_RC1.md`
Formal theorem certificate: `theorem_certificates/STAGE33_LEAN_FORMAL_VERIFICATION_CERTIFICATE.md`

## Substantive decisions explicitly attributable to the author (from conversation)

- Research question: public/private co-creation/innovation intermediation, distinguishing when a local government's free basic facility is economically warranted relative to private specialist matching.
- Research integrity: similar verbal conclusions in earlier literature do not by themselves disqualify a model with a different mechanism; nevertheless, closest papers must be cited and a distinct *economic* mechanism explained, not only a change of labels.
- Publication ambition: complete a bounded theoretical paper first, then determine suitable journal level from certified contribution rather than assume a high-tier target.
- Result selection: use existing Stage33 four-proposition structure and add lemmas/corollaries only where they clarify logically distinct economic steps.
- Formal method: explicitly requested Lean verification instead of relying solely on Z3, then instructed that Stage 8's canonical theory be finalized.

These are `AUTHOR` judgments about economic purpose, research direction, admissible literature positioning and theory-freeze intention. **They are not claimed to be personal proof verification.**

## Machine/AI verification actors

- `AI`: Stage 33 modelling support; Stage4 repair; numerical game analysis; Stage6 re-kill; Stage7 welfare accounting; Stage7.5A scope report; drafting formal source and frozen candidate.
- `COMPUTATION`: independent Python Cournot/payoff and price-deviation tests; 7x16 baseline audit, 13x16 welfare cases; exact rational coefficient regeneration on 2026-10-28.
- `FORMAL`: Lean 4 v4.32.1 and Mathlib at SHA `520045ab14e26149ee970e2e617ca04b09bde5d6`; CI success on exact formal source; 18 theorem statements accepted under encoded hypotheses, with proof-critical-core scope.
- `AUTHOR` proof verification: **NOT RECORDED**.
- `EXTERNAL HUMAN`: **NOT RECORDED**.

## Required author's own review — open, do not auto-sign

The author must personally inspect/correct the following economics/proof-level statements, and record substantive review tied to the exact frozen candidate:

1. **Price/sorting:** Why P-only `f=(U+T)/2` needs `U<T`; why full coverage uses `f=U`; why, with the other region fixed, `T` rising increases cutoff while total participant skill can fall. State what changes in P/P.
2. **Regional spillover:** Derive or critically inspect the foreign component `C ell[S_j^0+lambda K_j ell/4]`; interpret what negative `K_j` means for existing producers and the JV share; distinguish the two-profile comparison from the full planner optimum.
3. **Entry/game globality:** Explain why a zero of the B government addition gain only gives a conditional best-response switch; verify that the other three entrants' deviations and all non-equilibrium profiles must be examined and that at knife edge two pure NEs can coexist. Identify the price-uniqueness sufficient-condition domain and a regime outside it.
4. **Formal source fidelity:** Confirm that Lean proves the encoded polynomial signs, payoff identity and binary-flip logic from supplied primitive coefficients/assumptions, **not** the original Cournot derivation, full P/P game, global all16 SPNE or historical novelty.
5. **Publication bounds:** Accept the excluded claims (physical-facility necessity, subsidy dominance, unrestricted first best, arbitrary taste/skill correlation, universal entry reversals) as explicit limits.

A short author-written annotated proof sketch or documented critical correction is sufficient; the author need not rederive every Lean tactic or rewrite the entire mathematics. Chat instructions like `freeze it` are accepted as project authorization but **not** as evidence the author performed these five checks.

## Final status

`AUTHOR SIGN-OFF PENDING`. This record is a truthful basis for a later Stage 7.5A GO and Stage 8 Theory Freeze; it must not be modified to say PASS without an actual author review artifact. Formal-source commit and this exact candidate's SHA must be recorded at authentic sign-off time. No AI-generated signature is supplied.
