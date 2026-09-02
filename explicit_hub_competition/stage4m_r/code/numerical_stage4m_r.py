import numpy as np

N=500_000
rng=np.random.default_rng(20260902)
p=rng.uniform(0.01,0.99,N)
l=rng.uniform(0.02,1.0,N)
m=l+rng.uniform(0.02,1.5,N)
h=m+rng.uniform(0.02,1.5,N)
a=h+rng.uniform(0.02,5.0,N)
sigma=rng.uniform(0,5,N)
q=1-p

muS=p*(h+m)/2+p*q*l
nuS=p*(h*h+m*m)/2+p*q*l*l
muB=p*h+p*q*l+p*q*q*m
nuB=p*h*h+p*q*l*l+p*q*q*m*m
muB0=p*h+p*q*l
nuB0=p*h*h+p*q*l*l

def K(mu,nu):
    return 2*a*a/9+4*a*mu/9+11*nu/18-7*mu*mu/18

G=K(muB,nuB)-2*K(muS,nuS)+K(0,0)+sigma*(muB-2*muS)
G0=K(muB0,nuB0)-2*K(muS,nuS)+K(0,0)+sigma*(muB0-2*muS)
G_sigma0=K(muB,nuB)-2*K(muS,nuS)+K(0,0)
G_NS=(1+sigma)*(muB-2*muS)
G_NS0=(1+sigma)*(muB0-2*muS)
X1=(K(muB,nuB)-K(muS,nuS))+sigma*(muB-2*muS)

def counts(z):
    eps=1e-12
    return int((z>eps).sum()),int((z<-eps).sum()),int((np.abs(z)<=eps).sum())

print("referral Gamma +,-,0:",counts(G))
print("no-ref Gamma +,-,0:",counts(G0))
print("partner-surplus-zero Gamma +,-,0:",counts(G_sigma0))
print("no-downstream referral Gamma +,-,0:",counts(G_NS))
print("no-downstream no-ref Gamma +,-,0:",counts(G_NS0))
print("second-hub externality +,-,0:",counts(X1))
print("max Gamma:",float(G.max()))
print("max no-ref Gamma:",float(G0.max()))
print("max partner-zero Gamma:",float(G_sigma0.max()))
