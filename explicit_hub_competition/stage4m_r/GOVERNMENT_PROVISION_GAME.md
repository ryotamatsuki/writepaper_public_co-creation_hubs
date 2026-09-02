# Government Provision Game

Define moments for the single-hub distribution `S`:

`mu_S=p(h+m)/2 + p(1-p)l`

`nu_S=p(h^2+m^2)/2 + p(1-p)l^2`.

With both hubs and referral:

`mu_B=ph+p(1-p)l+p(1-p)^2m`

`nu_B=ph^2+p(1-p)l^2+p(1-p)^2m^2`.

Let `K_S=K(mu_S,nu_S)`, `K_B=K(mu_B,nu_B)`, `K_0=K(0,0)`.

For Region 1:

- `W_1(0,0)=K_0`
- `W_1(1,0)=K_S+2 sigma mu_S-F`
- `W_1(0,1)=K_S`
- `W_1(1,1)=K_B+sigma mu_B-F`.

The factor 2 in the single-own-hub partner term arises because the sole hub intermediates both regional projects; all matched partners are residents of its region.

Gross provision incentives are

`D_0^g=K_S-K_0+2 sigma mu_S`

`D_1^g=K_B-K_S+sigma mu_B`.

Net incentives subtract the same `F`.
