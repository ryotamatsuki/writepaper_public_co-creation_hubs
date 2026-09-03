"""Stage 7R-TP exact transfer-cancellation check."""
import sympy as sp
p,nT=sp.symbols("p nT", real=True)
regional_fee_term=-p*nT
private_profit=p*nT
assert sp.simplify(regional_fee_term+private_profit)==0
print("PASS: national private-fee transfer cancels exactly")
