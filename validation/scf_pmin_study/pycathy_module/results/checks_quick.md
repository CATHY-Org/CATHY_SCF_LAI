## quick

| check | result | detail |
|---|---|---|
| runs & water budget closes: q_scf_scf0_pmin-150 | PASS | completed=True residual=0.0092 % of V0 |
| runs & water budget closes: q_scf_scf0.3_pmin-150 | PASS | completed=True residual=0.0097 % of V0 |
| runs & water budget closes: q_scf_scf1_pmin-1e35 | PASS | completed=True residual=0.0094 % of V0 |
| runs & water budget closes: q_classic_scf0_pmin-150 | PASS | completed=True residual=0.0092 % of V0 |
| runs & water budget closes: q_classic_scf1_pmin-1e35 | PASS | completed=True residual=0.0094 % of V0 |
| SCF=0 -> no transpiration | PASS | max T = 0.00e+00 mm/d |
| SCF=1, PMIN off -> no evaporation | PASS | max E = 0.000 mm/d |
| SCF=0.3: E ~ 0.7*PET right after the pulse | PASS | E/PET = 0.700 |
| SCF=0.3: 0.15 PET < T <= 0.3 PET | PASS | T/PET = 0.281 |
| SCF=0.3: E + T <= PET | PASS | max (E+T)/PET = 0.984 |
| SCF build == classic at SCF=0 (cumulative E, 2 %) | PASS | scf 3.587 mm vs classic 3.587 mm |
| SCF build == classic at SCF=1 (cumulative T, 2 %) | PASS | scf 4.354 mm vs classic 4.354 mm |
| cumulative E non-decreasing: PMIN -5 -> -150 -> -1e4 | PASS | {-5.0: 2.97, -150.0: 3.587, -10000.0: 4.08} |
| surface head never below PMIN=-5 | PASS | min = -5.0000 m |
| surface head never below PMIN=-150 | PASS | min = -150.0000 m |
| SCF=1: T independent of PMIN (1-day, 5 %) | PASS | 4.354 vs 4.355 mm |
