"""Stage 3R-PMS diagnostic: four-way participation witness.

Purpose: show that rational non-use plus H1/H2/HT can coexist robustly.
This script is NOT a candidate matching model and intentionally uses deterministic
provider indices so that participation feasibility is separated from matching specificity.

Python >=3.10; no third-party dependencies.
"""


def shares(mu_t=0.65, c_pub=0.18, c_t=0.28, n=400):
    counts = [0, 0, 0, 0]
    for ik in range(n):
        k = ik / (n - 1)
        mu1 = 0.75 - 0.50 * k
        mu2 = 0.25 + 0.50 * k
        for iz in range(n):
            z = iz / (n - 1)
            values = [0.0, z * mu1 - c_pub, z * mu2 - c_pub, z * mu_t - c_t]
            h = max(range(4), key=lambda j: values[j])
            counts[h] += 1
    total = n * n
    return tuple(x / total for x in counts)


if __name__ == "__main__":
    print("baseline", shares())
    for dmu in (-0.02, 0.0, 0.02):
        for dc in (-0.01, 0.0, 0.01):
            s = shares(mu_t=0.65 + dmu, c_pub=0.18 + dc, n=200)
            print("perturb", dmu, dc, s, "all_positive", all(x > 0 for x in s))
