"""Stage 3R-TP symbolic precheck.
Python >=3.11; SymPy >=1.12.
No random draws are used.
This script verifies only algebraic cutoff and implicit-price-response identities.
"""
import sympy as sp

z, v, alpha, beta = sp.symbols('z v alpha beta', positive=True)
rho, rhoT, xi, Fi, FT = sp.symbols('rho rhoT xi Fi FT', real=True)
k_h, k_k, p_h, p_k = sp.symbols('k_h k_k p_h p_k', real=True)
Ph, Pk = sp.symbols('P_h P_k', positive=True)

Pi = sp.simplify(rho + xi + beta*Fi)
PT = sp.simplify(rhoT + beta*FT)

z_h0 = sp.solve(sp.Eq(z*(v + alpha*Ph) - k_h - p_h, 0), z)[0]
z_hk = sp.solve(sp.Eq(z*(v + alpha*Ph)-k_h-p_h,
                         z*(v + alpha*Pk)-k_k-p_k), z)[0]

D, Dp, Dpp, Dx, Dpx, p = sp.symbols('D Dp Dpp Dx Dpx p', real=True)
# From F(p,x)=D+p D_p=0, F_p=2D_p+pD_pp, F_x=D_x+pD_px.
p_x = sp.simplify(-(Dx+p*Dpx)/(2*Dp+p*Dpp))

print('P_i =', Pi)
print('P_T =', PT)
print('z_h0 =', sp.factor(z_h0))
print('z_hk =', sp.factor(z_hk))
print('p_x IFT =', p_x)
