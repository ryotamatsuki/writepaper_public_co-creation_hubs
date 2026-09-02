# Firm Choice Architecture

For a project of type `theta`, define ex ante expected net productive payoffs

- `U_1(theta)=M_1(theta,x_1)`;
- `U_2(theta)=M_2(theta,x_2)`;
- `U_T(theta)=M_T(theta)-p_T`;
- `U_0(theta)=0`.

The baseline local public access price is normalized to zero; real public cost enters government welfare, not project utility.

The project chooses

`h(theta) in argmax {U_1(theta), U_2(theta), U_T(theta), 0}`.

Ties are measure-zero under the selected continuous type/value specification or resolved by a fixed nonstrategic rule.

## Why this is not Hotelling

Heterogeneity comes from `v(theta,s)` and partner capability composition, not a taste-distance shock. Two projects can rank hubs differently because their productive requirements fit different capability sets.

## Why `0` is essential

Without the no-hub option, all projects are forced into an intermediary and extensive-margin additionality disappears. The baseline must allow a positive-measure set for which every available intermediary yields nonpositive expected net surplus.

## Ex ante choice, ex post matching

Projects choose using expected productive value. The realized partner and realized match surplus can differ. This separation is needed to define both intermediary choice and actual productive match creation without introducing an ad hoc taste term.