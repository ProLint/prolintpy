# Using `prolinpy` for protein-ligand interactions

The analysis and visualization workflow implemented by `prolintpy` can also be used to study protein-ligand interactions.

If the ligand is free to move then you can use `prolintpy` to visualize the different interactions that are formed during
the simulation. If the ligand is bound to the protein, you can still use the software to analyze protein-ligand interactions
by adjusting the distance cutoff. You will get a detailed profile of which residues are interacting with the ligand.

The 2D density app will probably not work for such usecases, and the 3D densities will be partially correct if the ligand is bound
to the protein, but these shortcomings can be corrected by centering the protein first.
