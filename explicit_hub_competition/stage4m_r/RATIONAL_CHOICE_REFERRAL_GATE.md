# Rational Choice × Referral Gate

## Comparison

- R1 full ex-ante compatibility information: rejected. Search/matching becomes largely mechanical and referral has no information role.
- **R2 project type known, partner-specific compatibility unknown: selected.**
- R3 broad specialization only: feasible but less minimal because it adds unnecessary uncertainty about project type/need.

## Selected finite structure

Each hub `i` has a resident generalist `g_i` and specialist `s_i`. Projects have type `theta in {1,2}`. For type 1:

- `v(theta1,s1)=h`
- `v(theta1,s2)=m`
- `v(theta1,g1)=v(theta1,g2)=l`

with `0<l<m<h`. Type 2 is symmetric. Each feasible partner is independently compatible with probability `p`.

The selected hub implements the highest-value compatible direct match. Only if all direct candidates fail can it refer the project to the rival hub's **non-overlapping specialist**.

## Fatal-gate calculation

For a type-1 project starting at Hub 1:

- `h` with probability `p`
- `l` with probability `(1-p)p`
- referral match `m` with probability `(1-p)^2 p`
- zero otherwise.

Starting at Hub 2 reverses the `h` and `m` positions:

- `m` with probability `p`
- `l` with probability `(1-p)p`
- referral match `h` with probability `(1-p)^2 p`
- zero otherwise.

Hence choosing the type-matched hub shifts mass

`p[1-(1-p)^2] = p^2(2-p) > 0`

from `m` to `h`, leaving the `l` and zero probabilities unchanged. The preferred first hub therefore **first-order stochastically dominates** the rival-first route.

This ranking holds for any strictly increasing beneficiary payoff in match productivity, including the Cournot profit used later.

## Gate verdict

**PASS — ENDOGENOUS HUB CHOICE IS RATIONAL AND NON-COSMETIC.**

Referral does not erase hub choice because the direct stage can preempt a higher-value specialist that would otherwise only be reached after failure. This ordering follows from the frozen direct-search-then-referral architecture, not from a taste shock, fee, distance or arbitrary home preference.
