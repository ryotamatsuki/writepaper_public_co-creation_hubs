# Four-Way Coexistence Test

A simple diagnostic provider-index witness (used only to test participation, not matching specificity) sets

`k,z ~ U[0,1]`,
`mu_1=0.75-0.50k`,
`mu_2=0.25+0.50k`,
`mu_T=0.65`,
public cost to the project `0.18`, metro project cost `0.28`.

Project choice maximizes `{0, z mu_1-0.18, z mu_2-0.18, z mu_T-0.28}`.

A dense-grid calculation gives shares near `(0.293,0.335,0.335,0.038)` for `(0,H1,H2,HT)`. Perturbing metro success by ±0.02 and public cost by ±0.01 in nine combinations leaves all four shares positive.

Therefore four-way coexistence has an open-set witness. The witness intentionally fails matching-specificity; its role is to isolate the problem. **Participation can be fixed without dynamics.**