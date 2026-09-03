# Participation Fixed Point

In the central regime define `delta=alpha beta`, `A_i=v+alpha(rho+x_i)`, `A_T=v+alpha rhoT`, `a=kappa_T+p_T`, `d=kappa_L-a`.

Then
`q_i=A_i+delta(1-t_i)`,
`q_T=A_T+delta(t_1+t_2-2a/q_T)`,
and the cutoffs solve
`t_i(q_i-q_T)=d`.

Thus the fixed point is the three-equation system
`F_i=t_i[A_i+delta(1-t_i)-q_T]-d=0`, i=1,2,
`F_T=q_T-A_T-delta[t_1+t_2-2a/q_T]=0`.

This is the minimal smooth system used for implicit differentiation and the numerical certificate.