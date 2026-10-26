"""Stage 5: off-path price-subgame counterexamples in the quadratic mixed Hotelling model.

This executable audit shows price-equilibrium multiplicity near collocation.
It does not prove either SPNE existence or SPNE nonexistence in the full game.
"""
from math import isclose

T, V, M = 1.0, 2.0, 1.0


def diagnose(delta, a, original_b, dev_b, price_G=1.0):
    assert 0 <= a < dev_b <= 1
    D = dev_b - a
    advantage_at_zero = delta - T * D * (a + dev_b)
    assert advantage_at_zero > 0
    price_P = price_G + advantage_at_zero
    assert price_P >= 2*T*D  # private profit maximized at all-P kink
    assert V-price_G-T*max(a*a, (1-a)**2) > 0  # outside option irrelevant

    candidate_profit = M*(T+4*delta)**2/(4*T)

    def demand_p(p):
        z = (a + dev_b)/2 + (p-price_G-delta)/(2*T*D)
        return min(1.0, max(0.0, 1.0-z))

    def profit_p(p):
        return M*p*demand_p(p) if p >= 0 else float("-inf")

    sampled_max = max((profit_p(i/10000), i/10000)
                      for i in range(30001))
    assert sampled_max[0] <= profit_p(price_P) + 1e-5
    # The public price is a best response: everyone is covered, and the
    # socially better private hub serves all types x>0.

    u_diff_at_zero = (delta-(price_P-price_G)
                      -T*((0-dev_b)**2-(0-a)**2))
    assert isclose(u_diff_at_zero, 0, abs_tol=1e-12)
    assert profit_p(price_P) > candidate_profit

    return {
        "delta": delta, "a": a,
        "b_original": original_b, "b_deviation": dev_b,
        "Delta": advantage_at_zero,
        "p_G": price_G, "p_P": price_P,
        "original_private_profit": candidate_profit,
        "deviation_private_profit": profit_p(price_P),
        "gain": profit_p(price_P)-candidate_profit,
        "max_safe_public_price": V-T*max(a*a,(1-a)**2),
        "sampled_private_best_reply": sampled_max,
    }


def low_continuation_profit(a, b, delta, t=1):
    """Private profit with the lowest feasible all-P price equilibrium.

    In an interior efficient split, the equilibrium is p_G=p_P.
    In all-P regions, take the minimum nonnegative public price
    q=max(0, 2*t*D-Delta); P gets q+Delta. Other pricing
    continuations also exist; this is an explicit selection rule.
    """
    D=abs(b-a)
    if D < 1e-12:
        return delta
    if b > a:
        z=(a+b)/2-delta/(2*t*D)
        if z <= 0:
            Delta=delta-t*D*(a+b)
            return max(Delta, 2*t*D)
        if z >= 1:
            return 0.0
        return 2*t*D*(1-z)**2
    z=(a+b)/2+delta/(2*t*D)
    if z >= 1:
        Delta=delta-t*D*(2-a-b)
        return max(Delta, 2*t*D)
    if z <= 0:
        return 0.0
    return 2*t*D*z*z


if __name__ == "__main__":
    for args in ((.09,.16,.66,.17),(.15,.10,.60,.11)):
        evidence=diagnose(*args)
        print(evidence)
        delta,a,b_star,_=args
        prices=[low_continuation_profit(a,b/10000,delta)
                for b in range(10001)]
        best=max(prices)
        at_candidate=low_continuation_profit(a,b_star,delta)
        assert abs(best-at_candidate)<1e-7
        print("LOW-PRICE CONTINUATION: max private profit",
              best, "candidate profit", at_candidate)
    print("PASS: off-path multiplicity, profitable high-price continuation,")
    print("and grid best-response test under selected low-price continuation.")
