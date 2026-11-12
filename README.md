# Boron Existence Form Calculator
Author: Tang Xinyu, Qinghai University. MIT license.

## Corrected model (2026-11-12)
Total concentration is mol B/L. All fractions are fractions of B atoms, not fractions of anion molecules. The previous app omitted polyborate B-number weights; use this corrected revision for recalculation.

Let x=[B(OH)3], h=10^(-pH). Molecular concentrations are:
|Species|Concentration|Default log10 K|
|---|---|---|
|B(OH)3|x|—|
|B(OH)4-|Ka*x/h|-9.2|
|B3O3(OH)4-|K1*x^3/h|-7.29|
|B5O6(OH)4-|K2*x^5/h|-6.77|
|B4O5(OH)4^2-|K3*x^4/h^2|-14.5|
|B3O3(OH)5^2-|K4*x^3/h^2|-16.3|

The corresponding formation reactions from boric acid are:
- B(OH)3 + H2O = B(OH)4- + H+
- 3B(OH)3 = B3O3(OH)4- + H+ + 2H2O
- 5B(OH)3 = B5O6(OH)4- + H+ + 5H2O
- 4B(OH)3 = B4O5(OH)4^2- + 2H+ + 3H2O
- 3B(OH)3 + H2O = B3O3(OH)5^2- + 2H+

Solve CT = x + c1 + 3c2 + 5c3 + 4c4 + 3c5 by bisection on [0,CT].
Report fi=ni*ci/CT. The six fractions sum to one.

These defaults reproduce the revised author model; original literature provenance and medium-specific calibration of each constant have not been independently verified here. Concentrations approximate activities, water activity is taken as unity, and h=10^(-pH) is used within that approximation. Ionic-strength corrections, ion pairing and temperature dependence are not implemented. Predictions at high ionic strength are conditional on these assumptions.

## Reproduce without a GUI
Install NumPy, then run:
```python
from boron_speciation_core import calculate_speciation, SPECIES
r = calculate_speciation(0.5, [10.0])
for name, fraction in zip(SPECIES, r["fractions"][:, 0]):
    print(name, 100*fraction)
print(r["mass_balance_residual"])
```
Expected percentages, in species order: 4.054662, 25.583190, 2.564043, 0.005816, 42.735502, 25.056786.
The file corrected_lookup.csv contains 726 concentration/pH conditions for comparison.

## Web interface
Install requirements.txt and run python app.py (Python 3.10 recommended for the existing dependency pins).
k1..k5 correspond to B3O3(OH)5^2-, B4O5(OH)4^2-, B5O6(OH)4-, B3O3(OH)4-, B(OH)4- respectively.
Y values are cumulative sums in that order; B(OH)3_fraction is exported separately. Integrals over pH are mathematical summaries, not physical species abundances.
The previously linked Hugging Face deployment is separate and has not been verified as synchronized with this corrected GitHub code.
