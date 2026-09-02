import sympy as sp
z,qh,qg,ah,ag=sp.symbols('z qh qg ah ag', nonzero=True)
print(sp.solve(sp.Eq(z*qh-ah,0),z)[0])
print(sp.solve(sp.Eq(z*qh-ah,z*qg-ag),z)[0])
