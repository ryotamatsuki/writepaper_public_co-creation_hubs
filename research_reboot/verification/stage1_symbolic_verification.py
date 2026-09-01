from sympy import symbols, simplify, factor, diff, solve

NA, NB = symbols("N_A N_B", positive=True)
m0, m1 = symbols("m0 m1", nonnegative=True)
beta, delta, vL, vX, eta, kappa, F, q0 = symbols(
    "beta delta v_L v_X eta kappa F q_0", positive=True
)
b = symbols("b", nonnegative=True)

PA = NA * (NA - 1) / 2
QAB = NA * NB
Lambda = (1 + beta) - beta * delta * (m1 + m0)

W00 = vL * PA * m0 + beta * (
    vL * PA * m0 * (1 - delta * m0) + vX * QAB * q0
)
W10 = vL * PA * m1 + beta * (
    vL * PA * m1 * (1 - delta * m1) + vX * QAB * q0
) - F

DeltaH_from_primitives = simplify(W10 - W00)
DeltaH_claim = vL * PA * (m1 - m0) * Lambda - F
assert simplify(DeltaH_from_primitives - DeltaH_claim) == 0

PiH = vL * PA * (m1 - m0) - F
DA = beta * vL * PA * (m1 - m0) * (1 - delta * (m1 + m0))
assert simplify(DeltaH_claim - (PiH + DA)) == 0

A = beta * eta * vX * QAB
pipeline_gain = A * b - kappa * b**2 / 2
bstar = solve(diff(pipeline_gain, b), b)[0]
assert simplify(bstar - A / kappa) == 0
assert simplify(diff(pipeline_gain, b, 2) + kappa) == 0

RAB = simplify(pipeline_gain.subs(b, bstar))
RAB_claim = A**2 / (2 * kappa)
assert simplify(RAB - RAB_claim) == 0

Delta_bundle = simplify(DeltaH_claim + RAB)
F_H_star = simplify(vL * PA * (m1 - m0) * Lambda)
F_Hb_star = simplify(F_H_star + RAB)
assert simplify((F_Hb_star - F_H_star) - RAB) == 0

derivatives = {
    "dDelta_dNA": factor(diff(DeltaH_claim, NA)),
    "dDelta_dF": factor(diff(DeltaH_claim, F)),
    "dDelta_ddelta": factor(diff(DeltaH_claim, delta)),
    "dDelta_dm1": factor(diff(DeltaH_claim, m1)),
    "dDelta_dm0": factor(diff(DeltaH_claim, m0)),
    "dbstar_dNA": factor(diff(bstar, NA)),
    "dbstar_dNB": factor(diff(bstar, NB)),
    "dbstar_deta": factor(diff(bstar, eta)),
    "dbstar_dkappa": factor(diff(bstar, kappa)),
    "dR_dNA": factor(diff(RAB, NA)),
    "dR_dNB": factor(diff(RAB, NB)),
}

print("STAGE1_SYMBOLIC_VERIFICATION_PASS")
print("Delta_H =", factor(DeltaH_claim))
print("Pi_H + D_A identity = PASS")
print("b* =", factor(bstar))
print("pipeline SOC =", factor(diff(pipeline_gain, b, 2)))
print("R_AB =", factor(RAB))
print("F_H+b* - F_H* =", factor(F_Hb_star - F_H_star))
for name, value in derivatives.items():
    print(name, "=", value)

print("\nInterpretive checks (not algebraic proofs):")
print("- q(b)=q0+eta*b*h makes policy effect of b zero when h=0 by definition.")
print("- d pipeline gain/db at b=0 = A > 0 under positive primitives, so b*>0 is built into the benchmark.")
print("- R_AB is additively separable from Delta_H; it does not reduce delta or alter local redundancy.")
print("- Region B contributes N_B only; it has no objective or strategic choice in the baseline.")
