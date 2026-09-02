import sympy as sp
t1,t2,qT,A1,A2,AT,d,delta,a=sp.symbols('t1 t2 qT A1 A2 AT d delta a')
q1=A1+delta*(1-t1); q2=A2+delta*(1-t2)
F=sp.Matrix([t1*(q1-qT)-d,t2*(q2-qT)-d,qT-AT-delta*(t1+t2-2*a/qT)])
print('F=',F)
print('J=',sp.factor(F.jacobian([t1,t2,qT])))
