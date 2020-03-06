# ToyFold-1D (in python)

Translated from Rhiju Das' Toyfold-1D in matlab, aiming to conserve code syntax / variable names / function names as much as possible. Docstrings mostly blatant plagiarizations!

## some implementation notes:

### Matlab dir -> py layout:

design -> went to design.Design or utils/design_utils.py

conformations -> conformation.Conformation or utils/conformation_utils.py

scoring -> inside conformation.Conformation or utils/energy_utils.py

util -> utils/conformation_utils.py

unittests -> tests (incomplete), not updated to class structure yet

### small things:

'partner'-type arrays are now zero-indexed, which changes unpaired state to be -1 (like Contrafold syntax)

Renamed input vec 'partner' to 'constraint' to distinguish it as an initial input from the p arrays.

