"""Standard-library exact mixed-copula half-plane moments."""
def moments(U,h,D,f,kappa):
    """E[(1,theta,z) 1{U+h theta+D z>=f}], exact for Fraction inputs.

    Integrate the independent square by splitting at its two edge crossings.
    This is equivalent to exact polynomial integration over a clipped polygon.
    Diagonal ties at a point have zero mass; full/empty corner ties are included.
    """
    zero=U*0;one=zero+1
    assert D>0 and h>=0 and 0<=kappa<=1
    cut=max(zero,min(one,(f-U)/(D+h)))
    diag=(one-cut,(one-cut*cut)/2,(one-cut*cut)/2)
    if h==0:
        t=max(zero,min(one,(f-U)/D))
        indep=(one-t,(one-t)/2,(one-t*t)/2)
    else:
        cuts=sorted(set([zero,one]+[x for x in ((f-U-D)/h,(f-U)/h) if 0<x<1]))
        vals=[zero,zero,zero]
        for l,r in zip(cuts,cuts[1:]):
            mid=(l+r)/2;t=(f-U-h*mid)/D
            if t>=1:continue
            if t<=0:
                vals[0]+=r-l;vals[1]+=(r*r-l*l)/2;vals[2]+=(r-l)/2
            else:
                A=(f-U)/D; B=h/D
                vals[0]+=(1-A)*(r-l)+B*(r*r-l*l)/2
                vals[1]+=(1-A)*(r*r-l*l)/2+B*(r**3-l**3)/3
                vals[2]+=((1-A*A)*(r-l)+A*B*(r*r-l*l)-B*B*(r**3-l**3)/3)/2
        indep=tuple(vals)
    return tuple(kappa*x+(1-kappa)*y for x,y in zip(diag,indep))
