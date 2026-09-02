import sympy as sp

A, beta, kappa, rho, tau = sp.symbols("A beta kappa rho tau", positive=True, real=True)
ki, kj = sp.symbols("ki kj", positive=True, real=True)

m = 2 + beta
n = beta - 1
b = 2*m*n

ai = 9*ki - 2*m**2
aj = 9*kj - 2*m**2
Delta = sp.expand(ai*aj - b**2)

xi = sp.factor(2*m*A*(aj+b)/Delta)
xj = sp.factor(2*m*A*(ai+b)/Delta)

qi = sp.factor((A + m*xi + n*xj)/3)
qj = sp.factor((A + m*xj + n*xi)/3)

qi_alt = sp.factor(3*A*ki*(aj+b)/Delta)
qj_alt = sp.factor(3*A*kj*(ai+b)/Delta)
assert sp.simplify(qi-qi_alt) == 0
assert sp.simplify(qj-qj_alt) == 0

pi_i = sp.factor(qi**2 - ki*xi**2/2)
pi_j = sp.factor(qj**2 - kj*xj**2/2)
pi_i_closed = sp.factor(A**2*ki*ai*(aj+b)**2/Delta**2)
pi_j_closed = sp.factor(A**2*kj*aj*(ai+b)**2/Delta**2)
assert sp.simplify(pi_i-pi_i_closed) == 0
assert sp.simplify(pi_j-pi_j_closed) == 0

Q = sp.factor(qi+qj)
CS = sp.factor(Q**2/2)
WiP = pi_i
WiW = sp.factor(pi_i + CS/2)
SW = sp.factor(pi_i + pi_j + CS)

H = kappa
L = kappa-rho
M = kappa-rho+tau

def subs_obj(expr, K_i, K_j):
    return sp.factor(expr.subs({ki: K_i, kj: K_j}))

profiles = {
    "00": (H,H),
    "10": (L,M),
    "01": (M,L),
    "11": (L,L),
}

W_P = {p: subs_obj(WiP,*K) for p,K in profiles.items()}
W_W = {p: subs_obj(WiW,*K) for p,K in profiles.items()}
SWp = {p: subs_obj(SW,*K) for p,K in profiles.items()}

D0P = sp.factor(W_P["10"]-W_P["00"])
D1P = sp.factor(W_P["11"]-W_P["01"])
GammaP = sp.factor(D1P-D0P)

D0W = sp.factor(W_W["10"]-W_W["00"])
D1W = sp.factor(W_W["11"]-W_W["01"])
GammaW = sp.factor(D1W-D0W)

S0 = sp.factor(SWp["10"]-SWp["00"])
S1 = sp.factor(SWp["11"]-SWp["01"])

x_i_sym, x_j_sym = sp.symbols("x_i x_j", nonnegative=True)
q_i_pre = sp.factor((A + m*x_i_sym + n*x_j_sym)/3)
assert sp.diff(q_i_pre, x_j_sym) == n/3

# Engagement best-response slope.
assert sp.simplify(sp.diff((2*m*A+b*x_j_sym)/ai, x_j_sym) - b/ai) == 0

# Exact beta=0 regional-welfare deletion proof.
ell, d, t, y = sp.symbols("ell d t y", positive=True, real=True)
GammaW_beta0 = sp.factor(
    GammaW.subs(beta,0).subs({
        kappa: ell+d+t,
        rho: d+t,
        tau: t
    })
)
GammaW_beta0_shift = sp.factor(GammaW_beta0.subs(ell, sp.Rational(4,3)+y))

positive = {
    A: 1, beta: sp.Rational(6,5), kappa: 5,
    rho: sp.Rational(3,2), tau: sp.Rational(6,5)
}
negative = {
    A: 1, beta: sp.Rational(1,2), kappa: 5,
    rho: sp.Rational(3,2), tau: sp.Rational(6,5)
}

print("q_i before engagement optimization:", q_i_pre)
print("dq_i/dx_j:", sp.diff(q_i_pre,x_j_sym))
print("engagement BR slope:", sp.factor(b/ai))
print("x_i*:", xi)
print("q_i*:", qi_alt)
print("pi_i*:", pi_i_closed)
print("GammaW beta=0 shifted:", GammaW_beta0_shift)
print("GammaW positive exact:", sp.factor(GammaW.subs(positive)))
print("GammaW negative exact:", sp.factor(GammaW.subs(negative)))
