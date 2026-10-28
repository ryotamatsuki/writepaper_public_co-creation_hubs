# Stage 10 — Reviewer verifiability and chain-of-proof index

Reference authority: `../STAGE8_CANONICAL_THEORY_FREEZE_V1.md`; Lean authority `../../formal/Stage33Formal.lean`; 18 proof-critical theorems only. Numerical authority: `verification/stage33_reproduce.py`. Neither formal source nor freeze modified.

| Claim | Model definitions | In-text proof bridge | Independent check and limitation |
|---|---|---|---|
| P1 full-participation dW/dlambda | Eq. (Wi), (KSi), full m=1 q=1/2 | derivative C K_i/4; sum omega=1 | frozen conditional proof; Cournot quantities and full matching held fixed; no Lean full market encoding |
| P2 posted-fee sorting | U_i, T_i, truncated cutoff; eq. (monopolyfee) | df= dT/2; dell=U dT/(2T²); m=1-ell, q=(1-ell²)/2 | strict interior fixed opposite region; fee kink outside interval |
| C1 fixed-mass selection | cutoff and uniform density | q-m/2=ell(1-ell)/2 | NOT reoptimised-equilibrium comparison |
| P3 public P→Both | H_i direct-service gain; differences in m and q | delta L0=C ell; delta L1=C ell²/4; foreign W_j incidence | restricted local/national pairwise comparison; not equivalent to 16-profile planner |
| P4 entry ownership | 4 independent binary entry flips, 16 profile W and pi | numeric witness including (B,P) vs (B,B), pooled equal-mean condition | Stage9 CI Python regression; negative Beta/Bertrand portability; no all-domain SPNE theorem |
| P/P price condition | simultaneous fees f_A,f_B and 2 coupled participation cutoffs | direct global best-response required, not isolated monopoly FOC | numerical baseline grid only plus upstream sufficient conditions; global arbitrary parameter domain UNSOLVED |
| Formal claim | 18 frozen statements | statement-to-scope map in theorem certificate | not a theorem of complete 4-player SPNE |

**Potential referee blockers to monitor:** (1) exact consistency of Cournot surplus formula and price-equilibrium branches in all 16 histories; (2) closed-form global P/P BR beyond tested baseline; (3) literature theorem-level absorption; (4) exact parameter certificate mapping from Cournot primitives to frozen baseline; (5) explicit source-based negative Beta/Bertrand sensitivity reproducer in Stage11.

**Editorial scope lock:** do not replace conditional theorem labels by global ones; do not equate positive expected venture profit with a positive regional welfare sign; do not subtract private fees twice; do not treat a statement as formally verified if only hand-derived.
