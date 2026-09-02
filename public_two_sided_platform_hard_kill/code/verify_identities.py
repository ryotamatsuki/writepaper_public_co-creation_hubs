import sympy as sp
A,B,k,kT,alpha,Delta=sp.symbols('A B k kT alpha Delta', positive=True)
p=((k-kT)*A-kT*B)/(2*(A+B))
assert sp.simplify(sp.diff(p,A)-k*B/(2*(A+B)**2))==0
Ax=-alpha/Delta**2
px=sp.factor(sp.diff(p,A)*Ax)
print('beta-zero p_x =',px)
assert px.is_negative is True
