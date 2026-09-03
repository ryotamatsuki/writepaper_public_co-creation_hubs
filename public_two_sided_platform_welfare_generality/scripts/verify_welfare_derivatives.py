"""Stage 7R-TP derivative identities. Run with numpy/scipy."""
from welfare_b3_vs_g import BASE,pstar,W,profit,solve_symmetric
def d1(f,a,b,h=1e-5): return (f(a+h,b)-f(a-h,b))/(2*h)
xg,pg,xb=solve_symmetric()
fw2=lambda a,b: W(2,a,b,pstar(a,b,BASE),BASE)
fpi=lambda a,b: profit(pstar(a,b,BASE),a,b,BASE)
print("dW_rival/dx =",d1(fw2,xg,xg))
print("dPi/dx       =",d1(fpi,xg,xg))
print("sum           =",d1(fw2,xg,xg)+d1(fpi,xg,xg))
