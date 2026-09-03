import sympy as sp
A,B,k,kT=sp.symbols('A B k kT', positive=True)
p=((k-kT)*A-kT*B)/(2*(A+B))
print('beta0 p* =',sp.factor(p))
print('dp/dA =',sp.factor(sp.diff(p,A)))
