# Minimality and Objective Audit

| Primitive | Classification | Deletion result |
|---|---|---|
| shared core | ESSENTIAL | no common network value to thicken |
| `beta>0` | ESSENTIAL FOR COMPLEMENTARITY | `Gamma^G<0` everywhere when `beta=0` |
| two gateways | ESSENTIAL | no bilateral cross-entry object |
| cross-region gateway usability | ESSENTIAL FOR OVERLAP | at `tau=rho` overlap vanishes and `Gamma^G>0` for `beta>0` |
| `tau` | ESSENTIAL FOR TRADE-OFF | `tau=0` makes gateways perfectly interchangeable and `Gamma^G<0` |
| `rho` | ESSENTIAL | no gateway access improvement without it |
| endogenous participation / outside option | ESSENTIAL | otherwise `N=2` and gateway entry cannot thicken the core |
| fixed cost `F` | EQUILIBRIUM ONLY | not needed for `Gamma^G`, needed for entry-game regimes |
| public/regional welfare objective | ESSENTIAL FOR SIGN-SWITCH THEOREM | natural private volume analogues have fixed signs |
| uniform `c` | TRACTABILITY | not mechanism-identified here |
| linear `beta N` | TRACTABILITY / UNTESTED ROBUSTNESS | no robustness extension authorized in Stage 4R |

## Public-objective relevance diagnostic

Baseline regional government values origin-region firm welfare, including firms' ability to use another region's gateway. That is what makes neighbor access both an outside option to own provision and a source of common-network value.

Three simple non-price diagnostic objectives do **not** reproduce the sign reversal:

1. Private gateway operator values all users routed through its facility. Then sole-gateway volume is `h+x` and two-gateway own volume is `s`; cross-effect is
   `s-(h+x)=-(B+rho-tau)/q<0`.
2. Gateway operator values only own-region users of its facility. Cross-effect is
   `s-h=beta tau/q>0`.
3. Objective is incremental own-region connected headcount. Cross-effect is
   `(s-x)-(h-r)=-(rho-tau)/q<0`.

Thus reasonable quantity-only private objectives generate fixed strategic signs. The **regional welfare objective is ESSENTIAL for the verified complements-to-substitutes reversal**, while the public label is not being used merely for welfare after the fact.

## Entry-game topology

If `D1<D0` (substitutes):

- `F<D1`: unique full entry;
- `D1<F<D0`: two asymmetric one-entry equilibria;
- `F>D0`: unique no entry.

If `D0<D1` (complements):

- `F<D0`: unique full entry;
- `D0<F<D1`: coordination multiplicity `(0,0)` and `(1,1)`;
- `F>D1`: unique no entry.