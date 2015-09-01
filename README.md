# wc-folder
A protein folding algorithm that maximizes entropy associated with water--e.g., minimizing hydrophobic surface area, maximizing hydrophilic surface area, and minimizing overall surface area.

Requires:

ECSPY (for PSO),
PyRosetta (if you want to use it for getting sequence from PDB, which you don't),
Biopython,
Pybel,
PeptideBuilder,
Numpy,
MSMS (executables)

Sorry.

-Your msms executables should be in a folder named msms, and you need to supply the path to that folder to the constructor.

e.g. if located at /home/msms give it /home.

-When you get the sequence another way, get rid of the code under "# Rosetta inits", except for self.sequence, which you want to be your sequence.

-Put Geometry.py from PeptideBuilder into the folder the pso/mc class is in.

See the example script for how to use the classes. Pretty self-explanatory.
