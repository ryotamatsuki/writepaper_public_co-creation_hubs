# Welfare Identity Check

Exact identity:

`W_1+W_2+Pi_T = [project gross/network benefits - real access costs] + [partner surplus] - [public investment costs]`.

Proof: the two regional project-surplus terms contain `-p_T[(t_1-s)+(t_2-s)] = -p_T n_T^F`; private profit contributes `+p_T n_T^F`.

The reproducibility script `scripts/verify_welfare_identities.py` checks the symbolic cancellation and the numerical equality at the frozen witness to machine precision.
