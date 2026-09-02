# Project Choice Partitions

Write gross intermediary values before the common real use cost as

`V_1(theta)=x_1[d-u theta]`,
`V_2(theta)=x_2[c+u theta]`,
`V_T=q:=m-p_T`.

A project chooses the largest of `V_1,V_2,q,kappa`; `kappa` is the gross-value threshold associated with the zero payoff from non-participation.

## Public-public cutoff

`theta_12=[x_1 d-x_2 c]/[u(x_1+x_2)]`.

When this lies in `(0,1)`, the public upper envelope is V-shaped and its minimum is

`K=a x_1 x_2/(x_1+x_2)`, `a=2c+u`.

## Public-metro cutoffs

For a constant metro gross offer `q`,

`theta_1T=[x_1 d-q]/(x_1 u)`,
`theta_2T=[q-x_2 c]/(x_2 u)`.

The exact clipped boundaries are

`L(q)=clip(theta_1T,0,1)`,
`R(q)=clip(theta_2T,0,1)`.

If `q>kappa`, metro demand occupies `[L(q),R(q)]` whenever `R(q)>L(q)`.

## Public-zero cutoffs

Replace `q` by `kappa` in the same formulas. When `q<=kappa`, metro cannot strictly beat non-participation; any middle gap is a zero-use interval.

## Key structural implication

A positive type-independent metro payoff (`q>kappa`) strictly dominates zero for every project. Therefore active metro use and a no-use interval cannot coexist in this baseline except at a measure-zero tie. This is why the active full game cannot simultaneously generate metro competition and extensive `0 -> H_i` creation.