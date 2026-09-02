import sympy as sp

x1, x2, c, u, m, p = sp.symbols('x1 x2 c u m p', positive=True)
A = (1/u) * (1/x1 + 1/x2)
B = 1 + 2*c/u
D = 2 * (A*(m-p) - B)
Pi = sp.expand(p*D)

p_star = sp.simplify(sp.solve(sp.diff(Pi, p), p)[0])
assert sp.simplify(p_star - (m/2 - (2*c+u)*x1*x2/(2*(x1+x2)))) == 0
assert sp.simplify(sp.diff(Pi, p, 2) + 4*A) == 0

dp_dx1 = sp.factor(sp.diff(p_star, x1))
dp_dx2 = sp.factor(sp.diff(p_star, x2))

theta_pp = sp.simplify((x1*(c+u)-x2*c)/(u*(x1+x2)))

print('D_T =', D)
print('p* =', p_star)
print('dp*/dx1 =', dp_dx1)
print('dp*/dx2 =', dp_dx2)
print('theta_PP =', theta_pp)
