from __future__ import annotations

from sympy import Rational, diff, expand, factor, simplify, sqrt, symbols


def check_zero(name: str, expr) -> None:
    reduced = simplify(factor(expand(expr)))
    if reduced != 0:
        raise AssertionError(f"{name} failed:\n{reduced}")
    print(f"[ok] {name}")


def check_positive_numeric(name: str, expr, samples: list[dict]) -> None:
    for index, subs in enumerate(samples, start=1):
        value = expr.subs(subs)
        value = simplify(value)
        if float(value) <= 0:
            raise AssertionError(
                f"{name} failed on sample {index}: value={value}, subs={subs}"
            )
    print(f"[ok] {name}")


def base_symbols():
    NA, NB = symbols("NA NB", positive=True)
    m0, m1 = symbols("m0 m1", positive=True)
    vL, vX = symbols("vL vX", positive=True)
    beta, delta, xi = symbols("beta delta xi", positive=True)
    eta0, eta1, q0 = symbols("eta0 eta1 q0", positive=True)
    F, kappa, b, h = symbols("F kappa b h", positive=True)
    return NA, NB, m0, m1, vL, vX, beta, delta, xi, eta0, eta1, q0, F, kappa, b, h


def base_samples():
    return [
        {
            "NA": Rational(6),
            "NB": Rational(3),
            "m0": Rational(1, 5),
            "m1": Rational(3, 5),
            "vL": Rational(2),
            "vX": Rational(3, 2),
            "beta": Rational(4, 5),
            "delta": Rational(3, 5),
            "xi": Rational(1, 10),
            "eta0": Rational(1, 20),
            "eta1": Rational(3, 20),
            "q0": Rational(1, 10),
            "F": Rational(5),
            "kappa": Rational(2),
        },
        {
            "NA": Rational(7),
            "NB": Rational(2),
            "m0": Rational(1, 10),
            "m1": Rational(11, 20),
            "vL": Rational(3, 2),
            "vX": Rational(5, 4),
            "beta": Rational(9, 10),
            "delta": Rational(2, 3),
            "xi": Rational(1, 8),
            "eta0": Rational(1, 25),
            "eta1": Rational(1, 10),
            "q0": Rational(2, 25),
            "F": Rational(4),
            "kappa": Rational(3),
        },
        {
            "NA": Rational(5),
            "NB": Rational(4),
            "m0": Rational(3, 20),
            "m1": Rational(1, 2),
            "vL": Rational(7, 5),
            "vX": Rational(6, 5),
            "beta": Rational(3, 4),
            "delta": Rational(7, 10),
            "xi": Rational(1, 12),
            "eta0": Rational(1, 15),
            "eta1": Rational(2, 15),
            "q0": Rational(1, 12),
            "F": Rational(7, 2),
            "kappa": Rational(5, 2),
        },
    ]


def sample_substitutions(symbol_map: dict[str, object]) -> list[dict]:
    substitutions = []
    for sample in base_samples():
        substitutions.append({symbol_map[name]: value for name, value in sample.items()})
    return substitutions


def verify_sections_03_to_06():
    (
        NA,
        NB,
        m0,
        m1,
        vL,
        vX,
        beta,
        delta,
        xi,
        eta0,
        eta1,
        q0,
        F,
        kappa,
        b,
        h,
    ) = base_symbols()

    PA = NA * (NA - 1) / 2
    QAB = NA * NB
    mh = m0 + (m1 - m0) * h
    Lambda = (1 + beta) - beta * delta * (m1 + m0)

    B1 = vL * PA * mh
    B2L = beta * vL * PA * mh * (1 - delta * mh + xi * b)
    B2X = beta * vX * QAB * (q0 + (eta0 + eta1 * h) * b)
    W = B1 + B2L + B2X - F * h - kappa * b**2 / 2
    Pi = vL * PA * mh + beta * vX * QAB * (q0 + (eta0 + eta1 * h) * b) - F * h - kappa * b**2 / 2

    A_P = vX * QAB * (eta0 + eta1 * h)
    A_W = xi * vL * PA * mh + A_P
    barPi = vL * PA * mh + beta * vX * QAB * q0 - F * h
    barW = vL * PA * mh + beta * (vL * PA * mh * (1 - delta * mh) + vX * QAB * q0) - F * h

    print("== Section 03 / model setup ==")
    check_zero("Pi compact form", Pi - (barPi + beta * A_P * b - kappa * b**2 / 2))
    check_zero("W compact form", W - (barW + beta * A_W * b - kappa * b**2 / 2))
    check_zero("dynamic wedge inside A_h", A_W - A_P - xi * vL * PA * mh)

    W00 = simplify(W.subs({h: 0, b: 0}))
    W10 = simplify(W.subs({h: 1, b: 0}))
    DeltaH = simplify(W10 - W00)
    DeltaH_expected = vL * PA * (m1 - m0) * Lambda - F

    print("== Section 04 / hub-alone efficiency ==")
    check_zero("Delta_H factorization", DeltaH - DeltaH_expected)
    FH_star = vL * PA * (m1 - m0) * Lambda
    check_zero("F_H^* representation", DeltaH - (FH_star - F))
    delta_star = (1 + beta) / (beta * (m1 + m0)) - F / (beta * vL * PA * (m1 - m0) * (m1 + m0))
    check_zero(
        "delta_H^* rearrangement",
        beta * (m1 + m0) * (delta_star - delta) - DeltaH / (vL * PA * (m1 - m0)),
    )
    check_zero(
        "d Delta_H / d N_A",
        diff(DeltaH_expected, NA) - vL * (2 * NA - 1) * (m1 - m0) * Lambda / 2,
    )
    check_zero("d Delta_H / d F", diff(DeltaH_expected, F) + 1)
    check_zero(
        "d Delta_H / d delta",
        diff(DeltaH_expected, delta) + beta * vL * PA * (m1 - m0) * (m1 + m0),
    )
    check_zero(
        "d Delta_H / d m1",
        diff(DeltaH_expected, m1) - vL * PA * ((1 + beta) - 2 * beta * delta * m1),
    )
    check_zero(
        "d Delta_H / d m0",
        diff(DeltaH_expected, m0) + vL * PA * ((1 + beta) - 2 * beta * delta * m0),
    )

    print("== Section 05 / decentralized support vs planner ==")
    AP0 = simplify(A_P.subs(h, 0))
    AP1 = simplify(A_P.subs(h, 1))
    AW0 = simplify(A_W.subs(h, 0))
    AW1 = simplify(A_W.subs(h, 1))

    bP0 = beta * AP0 / kappa
    bP1 = beta * AP1 / kappa
    bW0 = beta * AW0 / kappa
    bW1 = beta * AW1 / kappa

    check_zero("b_P^*(0) FOC", simplify(diff(Pi.subs(h, 0), b)).subs(b, bP0))
    check_zero("b_P^*(1) FOC", simplify(diff(Pi.subs(h, 1), b)).subs(b, bP1))
    check_zero("b_W^*(0) FOC", simplify(diff(W.subs(h, 0), b)).subs(b, bW0))
    check_zero("b_W^*(1) FOC", simplify(diff(W.subs(h, 1), b)).subs(b, bW1))
    check_zero("renewal gap h=0", bW0 - bP0 - beta * xi * vL * PA * m0 / kappa)
    check_zero("renewal gap h=1", bW1 - bP1 - beta * xi * vL * PA * m1 / kappa)
    check_zero(
        "renewal gap larger with hub",
        (bW1 - bP1) - (bW0 - bP0) - beta * xi * vL * PA * (m1 - m0) / kappa,
    )

    Pi0_star = simplify(Pi.subs({h: 0, b: bP0}))
    Pi1_star = simplify(Pi.subs({h: 1, b: bP1}))
    W0_star = simplify(W.subs({h: 0, b: bW0}))
    W1_star = simplify(W.subs({h: 1, b: bW1}))
    R0P = beta**2 * AP0**2 / (2 * kappa)
    R1P = beta**2 * AP1**2 / (2 * kappa)
    R0W = beta**2 * AW0**2 / (2 * kappa)
    R1W = beta**2 * AW1**2 / (2 * kappa)
    check_zero("Pi^*(0) quadratic value", Pi0_star - (barPi.subs(h, 0) + R0P))
    check_zero("Pi^*(1) quadratic value", Pi1_star - (barPi.subs(h, 1) + R1P))
    check_zero("W^*(0) quadratic value", W0_star - (barW.subs(h, 0) + R0W))
    check_zero("W^*(1) quadratic value", W1_star - (barW.subs(h, 1) + R1W))

    SA = vL * PA * (m1 - m0)
    DA = beta * vL * PA * (m1 - m0) * (1 - delta * (m1 + m0))
    MA = simplify((R1W - R0W) - (R1P - R0P))
    check_zero("Delta_P^*", (Pi1_star - Pi0_star) - (SA - F + (R1P - R0P)))
    check_zero("Delta_W^* wedge", (W1_star - W0_star) - ((Pi1_star - Pi0_star) + DA + MA))

    MA_expected = beta**2 * xi * vL * PA * (
        (m1 - m0) * (xi * vL * PA * (m1 + m0) + 2 * vX * QAB * eta0)
        + 2 * vX * QAB * eta1 * m1
    ) / (2 * kappa)
    check_zero("M_A factorization", MA - MA_expected)

    FP_star = SA + (R1P - R0P)
    FW_star = FP_star + DA + MA
    check_zero("F_P^* threshold", (Pi1_star - Pi0_star) - (FP_star - F))
    check_zero("F_W^* threshold", (W1_star - W0_star) - (FW_star - F))

    print("== Section 06 / pipeline-alone vs bundle ==")
    R0 = R0W
    R1 = R1W
    check_zero("pipeline-alone value", W0_star - W00 - R0)
    check_zero("bundle over pipeline-alone", W1_star - W0_star - (DeltaH + R1 - R0))
    Fdagger = FH_star + (R1 - R0)
    check_zero("F_dagger representation", (W1_star - W0_star) - (Fdagger - F))
    check_zero(
        "R1 - R0 positive factorization",
        (R1 - R0) - beta**2 * (AW1 - AW0) * (AW1 + AW0) / (2 * kappa),
    )

    samples = sample_substitutions(
        {
            "NA": NA,
            "NB": NB,
            "m0": m0,
            "m1": m1,
            "vL": vL,
            "vX": vX,
            "beta": beta,
            "delta": delta,
            "xi": xi,
            "eta0": eta0,
            "eta1": eta1,
            "q0": q0,
            "F": F,
            "kappa": kappa,
        }
    )

    check_positive_numeric("M_A positive on admissible samples", MA_expected, samples)
    check_positive_numeric("R1 - R0 positive on admissible samples", R1 - R0, samples)
    check_positive_numeric("R0 positive on admissible samples", R0, samples)

    return {
        "NA": NA,
        "NB": NB,
        "m0": m0,
        "m1": m1,
        "vL": vL,
        "vX": vX,
        "beta": beta,
        "delta": delta,
        "xi": xi,
        "eta0": eta0,
        "eta1": eta1,
        "q0": q0,
        "F": F,
        "kappa": kappa,
        "PA": PA,
        "QAB": QAB,
        "Lambda": Lambda,
        "FH_star": FH_star,
        "SA": SA,
        "DA": DA,
        "R0": R0,
        "R1": R1,
    }


def verify_section_07(context: dict):
    NA = context["NA"]
    NB = context["NB"]
    m0 = context["m0"]
    m1 = context["m1"]
    vL = context["vL"]
    vX = context["vX"]
    beta = context["beta"]
    delta = context["delta"]
    xi = context["xi"]
    eta0 = context["eta0"]
    eta1 = context["eta1"]
    F = context["F"]
    kappa = context["kappa"]
    PA = context["PA"]
    QAB = context["QAB"]
    Lambda = context["Lambda"]
    SA = context["SA"]
    DA = context["DA"]
    R0 = context["R0"]
    R1 = context["R1"]

    sigma, vLH, vLL, vXH, vXL, alpha = symbols(
        "sigma vLH vLL vXH vXL alpha", positive=True
    )
    Delta_m = m1 - m0

    print("== Section 07 / extensions ==")
    A0_alpha = alpha * xi * vL * PA * m0 + vX * QAB * eta0
    A1_alpha = alpha * xi * vL * PA * m1 + vX * QAB * (eta0 + eta1)
    bPalpha0 = beta * A0_alpha / kappa
    bPalpha1 = beta * A1_alpha / kappa
    bW0 = beta * (xi * vL * PA * m0 + vX * QAB * eta0) / kappa
    bW1 = beta * (xi * vL * PA * m1 + vX * QAB * (eta0 + eta1)) / kappa

    check_zero(
        "partial internalization gap h=0",
        bW0 - bPalpha0 - beta * (1 - alpha) * xi * vL * PA * m0 / kappa,
    )
    check_zero(
        "partial internalization gap h=1",
        bW1 - bPalpha1 - beta * (1 - alpha) * xi * vL * PA * m1 / kappa,
    )
    check_zero(
        "d b_{P,alpha}^*(0) / d alpha",
        diff(bPalpha0, alpha) - beta * xi * vL * PA * m0 / kappa,
    )
    check_zero(
        "d b_{P,alpha}^*(1) / d alpha",
        diff(bPalpha1, alpha) - beta * xi * vL * PA * m1 / kappa,
    )

    R0alpha = beta**2 * A0_alpha**2 / (2 * kappa)
    R1alpha = beta**2 * A1_alpha**2 / (2 * kappa)
    FPalpha = SA + alpha * DA + (R1alpha - R0alpha)
    FWdagger = SA + DA + (R1 - R0)
    XA = xi**2 * vL**2 * PA**2 * (m1**2 - m0**2)
    YA = 2 * xi * vL * PA * vX * QAB * (eta0 * (m1 - m0) + eta1 * m1)
    Omega_alpha = FWdagger - FPalpha

    check_zero("planner threshold equals bundle threshold", FWdagger - (context["FH_star"] + (R1 - R0)))
    check_zero(
        "Omega_alpha representation",
        Omega_alpha - (1 - alpha) * (DA + beta**2 * ((1 + alpha) * XA + YA) / (2 * kappa)),
    )
    check_zero("Omega_1 equals zero", simplify(Omega_alpha.subs(alpha, 1)))
    check_zero(
        "d Omega_alpha / d alpha",
        diff(Omega_alpha, alpha) + DA + beta**2 * (2 * alpha * XA + YA) / (2 * kappa),
    )

    vbarL = vLL * (1 - sigma) ** 2 + vLH * (1 - (1 - sigma) ** 2)
    vbarL_expected = vLL + (2 * sigma - sigma**2) * (vLH - vLL)
    check_zero("vbarL simplification", vbarL - vbarL_expected)
    check_zero(
        "vbarL derivative",
        diff(vbarL_expected, sigma) - 2 * (1 - sigma) * (vLH - vLL),
    )

    vbarX = vXL * (1 - sigma) + vXH * sigma
    vbarX_expected = vXL + sigma * (vXH - vXL)
    check_zero("vbarX simplification", vbarX - vbarX_expected)
    check_zero("vbarX derivative", diff(vbarX_expected, sigma) - (vXH - vXL))

    DeltaH_sigma = vbarL_expected * PA * Delta_m * Lambda - F
    check_zero(
        "d Delta_H(sigma) / d sigma",
        diff(DeltaH_sigma, sigma)
        - 2 * (1 - sigma) * (vLH - vLL) * PA * Delta_m * Lambda,
    )

    A1_sigma = xi * vbarL_expected * PA * m1 + vbarX_expected * QAB * (eta0 + eta1)
    b_sigma = beta * A1_sigma / kappa
    check_zero(
        "d b*(sigma) / d sigma",
        diff(b_sigma, sigma)
        - beta
        * (
            xi * 2 * (1 - sigma) * (vLH - vLL) * PA * m1
            + (vXH - vXL) * QAB * (eta0 + eta1)
        )
        / kappa,
    )

    R1_sigma = beta**2 * A1_sigma**2 / (2 * kappa)
    check_zero(
        "d R1(sigma) / d sigma",
        diff(R1_sigma, sigma)
        - beta**2
        * A1_sigma
        * (
            xi * 2 * (1 - sigma) * (vLH - vLL) * PA * m1
            + (vXH - vXL) * QAB * (eta0 + eta1)
        )
        / kappa,
    )


def verify_section_08_figure(context: dict):
    NA = context["NA"]
    NB = context["NB"]
    m0 = context["m0"]
    m1 = context["m1"]
    vL = context["vL"]
    vX = context["vX"]
    beta = context["beta"]
    delta = context["delta"]
    xi = context["xi"]
    eta0 = context["eta0"]
    eta1 = context["eta1"]
    kappa = context["kappa"]
    PA = context["PA"]
    QAB = context["QAB"]
    FH_star = context["FH_star"]

    print("== Section 08 / parameter-region figure ==")
    AW0 = xi * vL * PA * m0 + vX * QAB * eta0
    AW1 = xi * vL * PA * m1 + vX * QAB * (eta0 + eta1)
    R0 = beta**2 * AW0**2 / (2 * kappa)
    R1 = beta**2 * AW1**2 / (2 * kappa)
    Fdagger = FH_star + (R1 - R0)

    figure_formula = FH_star + beta**2 * (AW1**2 - AW0**2) / (2 * kappa)
    check_zero("Figure formula equals F_dagger", figure_formula - Fdagger)


def main():
    print("Python / SymPy verification for the Review of Regional Research working draft")
    context = verify_sections_03_to_06()
    verify_section_07(context)
    verify_section_08_figure(context)
    print("All symbolic and numeric checks passed.")


if __name__ == "__main__":
    main()
