# Journal of Economics: independent repricing and regular-distribution verification

Date: 2026-11-04 JST; Draft PR #73; article: Private Innovation Intermediation, Regional Incidence, and Public Access.

## Scoped new results
New manuscript source: manuscript/sections/04_general_extension.tex. New supplemental proof: manuscript/appendix/appendix_b_extension.tex. Reproducible exact checks: reviews/je_extension_cert.py.

**Corollary 1 (fixed fee).** In the original uniform model, for any lambda0 on the compact admissible interval and under Proposition 1's strict share band, fix the P-only posted fee at the optimum fP(lambda0). At lambda0 its cutoff derivative becomes -eta*ell/D<0. The EXACT local public cost-ceiling derivative is

    (-eta*ell/D)*[N_i*(u_i+d*ell)+C*(S_i^0+lambda0*K_i*ell/2)]
      +C*K_i*ell^2/4.

For sufficiently small positive gamma, the bracket is strictly positive, and K_i<0. Hence the fee-frozen slope is negative, while the equilibrium-fee slope is positive by Proposition 1. The price-freezing exercise is an out-of-equilibrium counterfactual, not a second provider Nash equilibrium. No parameters were retuned.

At original rational primitives, T_B(0)=0.040425071007120596; T_B(0.5)=0.04053167010244885 with optimal repricing, a +0.00010659909532825 increase. Freezing fee at the lambda=0 optimum gives T_B_frozen(0.5)=0.04000201059738977, a -0.00042306040973083 change.

**Proposition 2 (regular F).** At the zero-matching interior private fee optimum ell0, with local density g, positive density, nondecreasing hazard, and positive strict second-order denominator B_F, define J_F=integral[0,ell0] t*g(t)dt and A_F=u*(u+d*ell0)*g(ell0)^2/(d*B_F), where B_F=2*d*g(ell0)+(u+d*ell0)*g'(ell0). With full partner coverage and partner mean skill theta_bar_j,

    T_F'(lambda)=gamma*N_i*N_j*theta_bar_j*
      [M_i*J_F+omega_i*Pi*(J_F+A_F)]+O(gamma^2).

This has a positive weak-matching slope with K_i<0 and S_i^0<0 if
    (-M_i/Pi)*J_F/(J_F+A_F) < omega_i < min(1,-M_i/Pi)
is feasible. The direct-service gain must independently be positive for a positive gross cost-ceiling level. The uniform F, partner mean 1/2 reduction matches the original Proposition 1 IDENTICALLY in SymPy. No all-distribution or all-four-player equilibrium claim follows.

**Finite Beta(2,2).** At original gamma=0.1, omega_B=0.2, lambda in [0,0.5], region A fully participating, the Beta private-only price optimum is globally unique owing to increasing hazard. Exact Fraction FOC endpoint checks enclose its cutoff in (0.30,0.33) for every lambda in the interval. Monotonic rational lower bounds on the exact welfare derivative prove

    dT_B_Beta/dlambda > 0.0001372257882... > 1/10000

across the WHOLE closed lambda interval; this is not a two-point extrapolation. The previously reported failure of the FULL four-party Beta provider switch remains; no coupled off-path P/P equilibrium theorem is asserted for Beta.

## Predecessor scope
Damiano and Li (2007, 2008) have complementary fee selection; Kim, Li, Xu (2025) have participation distortions under entry fees; Ostrizek and Sartori (2023) have a general contract-screening/externalities model; Liu, Rivadeneyra, Reshidi (2026) establish private price responses to public provision. These ingredients remain pre-existing. Primary URLs respectively:
https://doi.org/10.1007/s00199-005-0058-2
https://doi.org/10.1162/JEEA.2008.6.4.789
https://doi.org/10.3982/TE5651
https://doi.org/10.1016/j.geb.2023.01.011
https://doi.org/10.34989/swp-2026-10

The existing proposition-level literature maps document failed direct restrictions, but Hua et al. technical results and certain final-version sources remain incomplete. This work does not certify general novelty or acceptance at Journal of Economics. Neither the new general proposition nor original Proposition 1 is part of the pinned limited Lean core. The complete four-party equilibrium baseline and frozen primitives are preserved.

## Status
These are newly proved conditional comparative statics in an explicitly authorised follow-up, not the previously frozen results. PR #73 remains Draft and unmerged. No submission is authorised until paper-specific author declarations and an active journal submission portal are confirmed.
