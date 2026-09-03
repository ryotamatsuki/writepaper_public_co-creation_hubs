# SOC and Denominator Audit

The implicit best-response formulas require negative own second derivatives. The frozen production result verifies

- `W_{i,x_ix_i}^{B3}<0` at the B3 witness;
- `W~_{i,x_ix_i}^{G}<0` at the G witness.

However, the historical counterexample search found nonconcave public stationary points in other primitive regions. Therefore neither SOC follows globally from the baseline primitive domain.

## Correct analytic domain

Define `Theta^{RI,SOC}` as the subset of the regular central branch satisfying:

1. participation Jacobian nonsingularity;
2. private FOC and `R_p<0`;
3. B3 public FOC and `W_{i,x_ix_i}^{B3}<0`;
4. G public FOC and `W~_{i,x_ix_i}^{G}<0`;
5. all strict interior/sorting inequalities.

Within this domain the sign of each best-response slope equals the sign of its corresponding cross derivative.

## Verdict

**Gate 3: PASS only by explicit SOC subregion, not by primitive-global proof.** This is legitimate for a local-branch theorem, but it materially weakens the prospect of a primitive Level-3A/B characterization because the SOC inequalities themselves remain equilibrium-dependent.