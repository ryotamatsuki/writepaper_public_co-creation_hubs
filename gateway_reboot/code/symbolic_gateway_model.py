"""Symbolic verification for Stage 3R-4R Gateway Reboot."""
import sympy as sp

B,beta,rho,tau,F = sp.symbols("B beta rho tau F", positive=True, real=True)
q = 1-2*beta

r = B/q
h = (B+rho-beta*tau)/q
x = h-tau
s = (B+rho)/q

D0 = sp.factor((h**2-r**2)/2)
D1 = sp.factor((s**2-x**2)/2)
Gamma = sp.factor(D1-D0)

X0 = sp.factor((x**2-r**2)/2)
X1 = sp.factor((s**2-h**2)/2)
G0 = sp.factor(D0+X0)
G1 = sp.factor(D1+X1)

Bstar = sp.factor(beta*(1-beta)*tau**2/(rho-tau) - (rho-tau)/2)

assert sp.simplify(
    Gamma - (
        2*beta*(1-beta)*tau**2
        -(rho-tau)*(2*B+rho-tau)
    )/(2*q**2)
) == 0
assert sp.simplify(Gamma.subs(B,Bstar)) == 0
assert sp.simplify((X1-X0)-Gamma) == 0
assert sp.simplify((G1-G0)-2*Gamma) == 0
assert sp.simplify((G0+G1)/2 - (D0+X1)) == 0
assert sp.simplify((G0+G1)/2 - (D1+X0)) == 0

print("r =", sp.factor(r))
print("x =", sp.factor(x))
print("h =", sp.factor(h))
print("s =", sp.factor(s))
print("D0 =", D0)
print("D1 =", D1)
print("Gamma =", Gamma)
print("X0 =", X0)
print("X1 =", X1)
print("G0 =", G0)
print("G1 =", G1)
print("B* =", Bstar)
print("dB*/dbeta =", sp.factor(sp.diff(Bstar,beta)))
print("dB*/drho =", sp.factor(sp.diff(Bstar,rho)))
print("dB*/dtau =", sp.factor(sp.diff(Bstar,tau)))

# Private-objective diagnostics (not baseline theory)
Q10_all = h+x
Q11_all = s
print("private total gateway-users cross-effect =", sp.factor(Q11_all-Q10_all))
print("private own-region gateway-users cross-effect =", sp.factor(s-h))
print("local connected-participation cross-effect =", sp.factor((s-x)-(h-r)))
