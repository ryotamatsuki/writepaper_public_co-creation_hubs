"""Stage 7R-TP phase-map logic."""
regions=[('fixed-price substitutes','M_B3 < 0'),('no fixed-price interaction','M_B3 = 0'),('complements survive pricing','M_B3 > 0 and M_G > 0'),('P2-R4','M_B3 > 0 > M_G')]
for r,c in regions: print(r,':',c)
print('Stage 4: beta=0 collapses relevant B3 cross effect; beta=.05 witness is P2-R4.')
