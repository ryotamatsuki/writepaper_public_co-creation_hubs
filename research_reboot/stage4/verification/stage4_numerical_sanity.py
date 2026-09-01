import math


def solve(v, B, A):
    s = math.sqrt(1 + 4*v)
    n0 = (s - 1)/2
    nR = (v*B - A)/(2*A)
    nSP = (v*B - A)/(2*A - v)
    alphaR = 1 - nR*(1+nR)/v
    alphaSP = 1 - nSP*(1+nSP)/v

    def n_of_alpha(alpha):
        return (math.sqrt(1 + 4*v*(1-alpha)) - 1)/2

    def UR(alpha):
        n = n_of_alpha(alpha)
        return A*alpha + B*n

    def W(alpha):
        n = n_of_alpha(alpha)
        return UR(alpha) + n*n/2

    return {
        'n0': n0,
        'nR': nR,
        'nSP': nSP,
        'alphaR': alphaR,
        'alphaSP': alphaSP,
        'UR_alphaR': UR(alphaR),
        'UR_1': UR(1.0),
        'W_alphaSP': W(alphaSP),
        'W_alphaR': W(alphaR),
    }


# Named benchmark.
benchmark = solve(v=1.0, B=1.0, A=0.8)
print('benchmark:', benchmark)
assert 0 < benchmark['alphaSP'] < benchmark['alphaR'] < 1
assert 0 < benchmark['nR'] < benchmark['nSP'] < 1
assert benchmark['UR_alphaR'] > benchmark['UR_1']
assert benchmark['W_alphaSP'] > benchmark['W_alphaR']

# Grid through an open collection of admissible parameter values.
count = 0
for v in [0.2, 0.5, 1.0, 1.5, 1.9]:
    for B in [0.6, 0.8, 1.0, 1.5, 2.0]:
        s = math.sqrt(1 + 4*v)
        lower = v/2 + v*(B-0.5)/s
        upper = v*B
        if not lower < upper:
            continue
        A = (lower + upper)/2
        out = solve(v=v, B=B, A=A)
        assert 0 < out['alphaSP'] < out['alphaR'] < 1
        assert 0 < out['nR'] < out['nSP'] < 1
        assert out['UR_alphaR'] > out['UR_1']
        assert out['W_alphaSP'] > out['W_alphaR']
        count += 1

print('PASS: admissible grid points =', count)
