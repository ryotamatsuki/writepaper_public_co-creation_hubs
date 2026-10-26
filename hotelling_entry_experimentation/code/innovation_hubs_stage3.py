"""Stage 3 reproducibility: binary pilot scope, endogenous investment, private entry.

This is a fixed-location (G=0, P=1) benchmark, NOT a solution of the
endogenous-location or arbitrary off-path pricing model. Units are normalized.
No third-party packages are required.
"""
from math import sqrt

T, V, ALPHA, SIGMA, GAMMA, B, K = 1.0, 2.0, .30, 1.40, 1.0, .50, .04
F_PRIVATE, F_PUBLIC = 1.0, .60


def pilot(n):
    r = ALPHA * n / (1 + SIGMA * (n - 1))
    investment = r  # argmax_e r e - e^2/2
    delta = GAMMA * n * investment
    value = r * investment - investment**2 / 2 + B * n * investment - K * n
    return investment, delta, value


def private_monopoly(mass):
    assert 0 < V < 3 * T
    q = sqrt(V / (3 * T))
    gross_profit = mass * (2 * V / 3) * q
    entered = gross_profit >= F_PRIVATE
    gross_surplus = mass * (V * q - T * q**3 / 3)
    return entered, gross_profit, gross_surplus - F_PRIVATE if entered else 0.0


def public_policy(mass, n):
    invest, delta, pilot_value = pilot(n)
    assert 0 <= delta < T
    gross_private_profit = mass * (T + delta)**2 / (2 * T)
    entry_surplus = mass * (T + delta)**2 / (4 * T)
    entered = gross_private_profit >= F_PRIVATE
    public_only_surplus = mass * (V - T / 3)
    welfare_before_public_fixed = public_only_surplus + pilot_value
    if entered:
        welfare_before_public_fixed += entry_surplus - F_PRIVATE
    return {
        "n": n, "investment": invest, "spillover": delta,
        "pilot_social_surplus": pilot_value,
        "private_entry": entered,
        "private_gross_profit": gross_private_profit,
        "entry_gross_social_surplus": entry_surplus,
        "welfare_before_public_fixed": welfare_before_public_fixed,
        "welfare_net": welfare_before_public_fixed - F_PUBLIC,
    }


def evaluate(mass):
    no_public = private_monopoly(mass)
    policies = [public_policy(mass, n) for n in (1, 2)]
    best = max(policies, key=lambda p: p['welfare_net'])
    build = best['welfare_net'] >= no_public[2]
    return best['n'], build, no_public, policies


if __name__ == '__main__':
    for mass in (0.30, 0.88, 0.90, 1.00, 1.18, 1.19, 1.50):
        scope, build, baseline, policies = evaluate(mass)
        print(f'M={mass:.2f} n*={scope} build={build} '
              f'no_public_W={baseline[2]:.6f} '
              + ' '.join(f'W{p["n"]}={p["welfare_net"]:.6f}' for p in policies))
    at_one = evaluate(1.0)
    assert at_one[0] == 1 and at_one[1]
    p1, p2 = at_one[3]
    assert not p1['private_entry'] and p2['private_entry']
    assert p2['pilot_social_surplus'] > p1['pilot_social_surplus']
    assert p2['entry_gross_social_surplus'] < F_PRIVATE
    print('thresholds: M2=', F_PRIVATE / ((T + pilot(2)[1])**2 / (2*T)),
          'M1=', F_PRIVATE / ((T + pilot(1)[1])**2 / (2*T)))
