import sympy as sp

alpha, A, B, v, nbar = sp.symbols(
    'alpha A B v nbar', positive=True, real=True
)

D = 1 + 4*v*(1-alpha)
n = (sp.sqrt(D)-1)/2

# Participation identities.
assert sp.simplify(n*(1+n) - v*(1-alpha)) == 0

dn = sp.simplify(sp.diff(n, alpha))
d2n = sp.simplify(sp.diff(n, alpha, 2))

assert sp.simplify(dn + v/sp.sqrt(D)) == 0
assert sp.simplify(d2n + 2*v**2/D**sp.Rational(3, 2)) == 0

# Regional sponsor.
U = sp.simplify(A*alpha + B*n)
dU = sp.simplify(sp.diff(U, alpha))
d2U = sp.simplify(sp.diff(U, alpha, 2))

assert sp.simplify(dU - (A - B*v/sp.sqrt(D))) == 0
assert sp.simplify(d2U + 2*B*v**2/D**sp.Rational(3, 2)) == 0

alpha_R = sp.simplify(1 + 1/(4*v) - B**2*v/(4*A**2))
assert sp.simplify(dU.subs(alpha, alpha_R)) == 0

n_R_formula = sp.simplify((v*B-A)/(2*A))
# Rather than simplify through square-root absolute values, verify the
# participation identity at the candidate n_R and alpha_R.
assert sp.simplify(n_R_formula*(1+n_R_formula) - v*(1-alpha_R)) == 0

# Planner.
S_X = sp.simplify(n**2/2)
W = sp.simplify(U + S_X)
dW = sp.simplify(sp.diff(W, alpha))
d2W = sp.simplify(sp.diff(W, alpha, 2))

planner_foc_target = A - v/sp.Integer(2) - v*(B-sp.Rational(1,2))/sp.sqrt(D)
assert sp.simplify(dW - planner_foc_target) == 0
assert sp.simplify(d2W - v**2*(1-2*B)/D**sp.Rational(3,2)) == 0
assert sp.simplify(dW - (dU + n*dn)) == 0

alpha_SP = sp.simplify(
    1 + 1/(4*v)
    - v*(B-sp.Rational(1,2))**2/(4*(A-v/sp.Integer(2))**2)
)

# Verify planner FOC by substituting the implied positive sqrt(D) relation,
# valid on the target region A>v/2 and B>1/2.
sqrtD_SP = v*(B-sp.Rational(1,2))/(A-v/sp.Integer(2))
assert sp.simplify(
    A-v/sp.Integer(2)-v*(B-sp.Rational(1,2))/sqrtD_SP
) == 0

n_SP_formula = sp.simplify((v*B-A)/(2*A-v))
assert sp.simplify(
    n_SP_formula*(1+n_SP_formula) - v*(1-alpha_SP)
) == 0

# Fixed-participation counterfactual.
U_F = sp.simplify(A*alpha + B*nbar)
outside_surplus_F = sp.simplify(
    nbar*v*(1-alpha)/(1+nbar) - nbar**2/2
)
W_F = sp.simplify(U_F + outside_surplus_F)

assert sp.simplify(sp.diff(U_F, alpha) - A) == 0
assert sp.simplify(
    sp.diff(W_F, alpha) - (A - nbar*v/(1+nbar))
) == 0

print('PASS: participation identity')
print('n(alpha) =', n)
print('dn/dalpha =', dn)
print('d2n/dalpha2 =', d2n)
print('PASS: regional FOC/SOC identities')
print('alpha_R =', alpha_R)
print('n_R =', n_R_formula)
print('PASS: planner FOC/SOC identities')
print('alpha_SP =', alpha_SP)
print('n_SP =', n_SP_formula)
print('PASS: dW = dU + n*dn decomposition')
print('PASS: fixed-participation counterfactual')
