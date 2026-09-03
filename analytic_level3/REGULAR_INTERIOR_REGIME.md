# Regular-Interior Regime

Let `delta=alpha beta`, `a=kappa_T+p`, `d=kappa_L-a`, `A_i=v+alpha(rho+x_i)`, and `A_T=v+alpha rho_T`.

The analytic track defines `Theta^RI` as an equilibrium branch satisfying all of the following. These conditions are deliberately separated by logical status.

## 1. Primitive restrictions

`v>0`, `alpha>0`, `beta>=0`, `rho>=0`, `rho_T>=0`, `gamma>0`, `kappa_L>0`, `kappa_T>=0`, `tau>=0`, and controls `x_i in (0,1)`, `p>0` for an interior private optimum.

## 2. Central sorting/interiority inequalities

For each home region `i`:

- `d=kappa_L-kappa_T-p>0`;
- `q_i>q_T>0`;
- `s=a/q_T` and `t_i=d/(q_i-q_T)` satisfy `0<s<t_i<1`;
- project masses `n_i^F=1-t_i>0` and `n_T^F=t_1+t_2-2s>0`;
- partner masses `b_i=rho+x_i+beta(1-t_i)` and `b_T=rho_T+beta n_T^F` lie in `(0,1)`;
- to retain the local home-public ordering against the rival public route, impose the corresponding strict no-rival-switch inequalities. At a symmetric point with `tau>0` these hold locally; off symmetry they must be checked rather than inferred from symmetry.

Strict inequalities exclude indifference masses and active-set boundaries.

## 3. Participation regularity

The fixed-point system in `(t_1,t_2,q_T)` must have a nonsingular Jacobian. With

`J11=q_1-q_T-delta t_1`, `J22=q_2-q_T-delta t_2`, `J13=-t_1`, `J23=-t_2`, `J31=J32=-delta`, `J33=1-2 delta a/q_T^2`,

require `det J != 0`.

The explicit quadratic reduction also requires a real positive branch:

`Delta_Q=[A_T+delta(t_1+t_2)]^2-8 delta a > 0` when `delta>0`, together with selection of the root that satisfies all central-regime inequalities. The algebraic lower root is not silently discarded; it is excluded only when it violates the stated regime conditions.

## 4. Private-price regularity and SOC

For `D_T=n_T^F` and `Pi_T=p D_T`, require the interior FOC

`R=D_T+p D_{T,p}=0`

and strict SOC/IFT denominator

`R_p=2D_{T,p}+pD_{T,pp}<0`.

This gives a locally unique smooth private-price policy conditional on the participation branch. It is not a global uniqueness claim.

## 5. Public FOCs, SOCs, and BR regularity

For each regime `m in {B3,G}`, require an interior public stationary point `F_i^m=0` and strict own curvature

`F_{i,x_i}^m<0`.

These SOCs are not consequences of the baseline primitive domain: the earlier counterexample search contains nonconcave public stationary points. Therefore the Level-3 domain is properly `Theta^{RI,SOC}`, the subset of the central regular branch satisfying both public SOCs and the required Jacobian nonsingularity. This condition is explicit rather than hidden.

## 6. Scope of uniqueness

The regime guarantees only local uniqueness of the selected participation branch, the interior private optimum, and each differentiable public best-response branch where the stated Jacobians/SOCs hold. Global equilibrium uniqueness is not assumed or claimed.