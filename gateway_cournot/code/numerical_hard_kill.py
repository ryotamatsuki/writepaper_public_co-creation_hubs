import numpy as np

SEED = 20260902
N = 500_000

def threshold(beta):
    m = 2.0 + beta
    n = beta - 1.0
    return (2*m*m + np.abs(2*m*n))/9.0

def profile_vec(A,beta,ki,kj):
    m=2.0+beta
    n=beta-1.0
    b=2.0*m*n
    ai=9.0*ki-2.0*m*m
    aj=9.0*kj-2.0*m*m
    Delta=ai*aj-b*b
    xi=2.0*m*A*(aj+b)/Delta
    xj=2.0*m*A*(ai+b)/Delta
    qi=3.0*A*ki*(aj+b)/Delta
    qj=3.0*A*kj*(ai+b)/Delta
    pi=A*A*ki*ai*(aj+b)**2/(Delta*Delta)
    pj=A*A*kj*aj*(ai+b)**2/(Delta*Delta)
    Q=qi+qj
    CS=0.5*Q*Q
    return {
        "xi":xi,"xj":xj,"qi":qi,"qj":qj,
        "pi":pi,"pj":pj,"Q":Q,"CS":CS,
        "WiP":pi,"WiW":pi+0.5*CS,
        "SW":pi+pj+CS,"Delta":Delta,"ai":ai,"aj":aj,
    }

rng=np.random.default_rng(SEED)
beta=rng.uniform(0.0,2.0,N)
A=rng.uniform(0.1,3.0,N)
L=threshold(beta)+rng.uniform(0.25,5.0,N)
rho=rng.uniform(0.05,2.5,N)
tau=rho*rng.uniform(0.02,0.98,N)
kappa=L+rho

H=kappa
M=L+tau
o00=profile_vec(A,beta,H,H)
o10=profile_vec(A,beta,L,M)
o01=profile_vec(A,beta,M,L)
o11=profile_vec(A,beta,L,L)

D0P=o10["WiP"]-o00["WiP"]
D1P=o11["WiP"]-o01["WiP"]
GammaP=D1P-D0P

D0W=o10["WiW"]-o00["WiW"]
D1W=o11["WiW"]-o01["WiW"]
GammaW=D1W-D0W

S0=o10["SW"]-o00["SW"]
S1=o11["SW"]-o01["SW"]

tol=1e-10
print("admissible draws:",N)
print("producer Gamma>0:",int(np.sum(GammaP>tol)))
print("producer Gamma<0:",int(np.sum(GammaP<-tol)))
print("producer near zero:",int(np.sum(np.abs(GammaP)<=tol)))
print("regional Gamma>0:",int(np.sum(GammaW>tol)))
print("regional Gamma<0:",int(np.sum(GammaW<-tol)))
print("regional near zero:",int(np.sum(np.abs(GammaW)<=tol)))

# Fixed-cost welfare-preview stress test.
rngF=np.random.default_rng(SEED+1)
gross_scale=np.maximum.reduce([
    D0W,D1W,np.abs(S0),np.abs(S1),np.full(N,1e-8)
])
F=rngF.uniform(0.0,1.25,N)*gross_scale

sub=D1W<D0W
comp=~sub
ne0=(sub&(F>D0W)) | (comp&(F>D0W))
ne2=(sub&(F<D1W)) | (comp&(F<D1W))
ne1=sub&(F>D1W)&(F<D0W)

V0=np.zeros(N)
V1=S0-F
V2=S0+S1-2*F
sp=np.argmax(np.vstack([V0,V1,V2]),axis=0)

print("no-entry NE / social one gateway:",int(np.sum(ne0&(sp==1))))
print("no-entry NE / social two gateways:",int(np.sum(ne0&(sp==2))))
print("full-entry NE / social one gateway:",int(np.sum(ne2&(sp==1))))
print("full-entry NE / social zero gateways:",int(np.sum(ne2&(sp==0))))
print("one-entry NE / social two gateways:",int(np.sum(ne1&(sp==2))))
print("one-entry NE / social zero gateways:",int(np.sum(ne1&(sp==0))))

# Independent direct linear-system checks for 10,000 draws.
check_rng=np.random.default_rng(1234)
idx=check_rng.choice(N,size=10_000,replace=False)
max_x_error=0.0
max_foc_residual=0.0

for z in idx:
    bz=float(beta[z]); Az=float(A[z]); kz=float(kappa[z])
    rz=float(rho[z]); tz=float(tau[z])
    mz=2.0+bz
    nz=bz-1.0
    bc=2.0*mz*nz
    Lz=kz-rz
    Mz=Lz+tz
    pairs=((kz,kz),(Lz,Mz),(Mz,Lz),(Lz,Lz))
    closed=(o00,o10,o01,o11)
    for p,(ki,kj) in enumerate(pairs):
        ai=9.0*ki-2.0*mz*mz
        aj=9.0*kj-2.0*mz*mz
        mat=np.array([[ai,-bc],[-bc,aj]])
        rhs=np.array([2.0*mz*Az,2.0*mz*Az])
        sol=np.linalg.solve(mat,rhs)
        x_closed=np.array([closed[p]["xi"][z],closed[p]["xj"][z]])
        max_x_error=max(max_x_error,float(np.max(np.abs(sol-x_closed))))
        max_foc_residual=max(
            max_foc_residual,
            float(np.max(np.abs(mat@x_closed-rhs)))
        )

print("max independent x error:",max_x_error)
print("max engagement FOC residual:",max_foc_residual)
