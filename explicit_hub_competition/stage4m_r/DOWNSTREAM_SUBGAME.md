# Downstream Cournot Subgame

Inverse demand is `P=A-Q`; baseline marginal cost is `c`; define `a=A-c`.

If beneficiary productivity gains are `x_i,x_j`, equilibrium quantities are

`q_i=(a+2x_i-x_j)/3`,
`q_j=(a+2x_j-x_i)/3`.

Profits are `pi_i=q_i^2`, and consumer surplus is

`CS=(q_i+q_j)^2/2`.

Regional market welfare uses local firm profit plus one half of common consumer surplus:

`w_i^M=pi_i+CS/2`.

For an iid productivity distribution with mean `mu` and second moment `nu`,

`K(mu,nu)=E[w_i^M]`

equals

`2a^2/9 + 4a mu/9 + 11 nu/18 - 7 mu^2/18`.

This single formula exactly covers all states `x_i,x_j in {0,l,m,h}`.

State enumeration over 10,000 random parameter vectors matches the moment formula with maximum absolute error `1.4210854715202004e-14`.
