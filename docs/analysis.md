# Analysis Reference

`prolintpy` automates analysis and visualization of lipid-protein interactions. <br>
This section outlines how to use ProLint for generating topologies and analyzing contact information.

### Loading Data Files
Trajectory IO processes are delegated to MDTraj and both
<a href="https://mdtraj.org/1.9.4/api/generated/mdtraj.Trajectory.html">MDTraj.Trajectory</a> and
<a href="https://mdtraj.org/1.9.3/api/generated/mdtraj.Topology.html">MDTraj.Topology</a>
objects are used to process the input files. `prolintpy` has its own topology classes that
are used to define both Proteins and Lipids in the system.

To load simulation files, you can use the following example code:

```python
import mdtraj as md
import prolintpy as pl

t = md.load('system.xtc', top='system.gro')
```

`prolintpy` assumes that the provided input files contain only proteins and lipids (i.e. no water or
ions particles). If ligands are not removed then they will be treated as lipids
(useful if you want to study protein-ligand interactions). If you require any trajectory manipulation (e.g. use a stride, remove
periodicity, concatenate trajectories) then you should do that before loading the files.

The reason why `prolintpy` does not do any manipulation of input trajectories is because there
is no reason to support such use-cases when they can be done using other tools.


### Topology definitions

Files loaded using MDTraj lack any topological description. It is easy to get lipids since they have a clear residue name in the
input coordinate files, but this is not the case for proteins. MD Simulation systems quite often contain complex bilayer systems with
different number of proteins embedded in them. In order to automate the analysis of lipid-protein interactions we have to build
a topological description of the system.

We start by defining the protein and lipid topologies:

```python
# Specify the resolution of the input data: martini or atomistic
resolution = "martini"

# Define the lipid topology
lipids = pl.Lipids(t.topology, resolution=resolution)

# Define the protein topology
proteins = pl.Proteins(t.topology, resolution=resolution).system_proteins(merge=True)
```

The method `system_proteins` will build the protein topology from a coordinate file. The `merge` option indicates if similar proteins should
be grouped or not. For instance, if a system contains four copies of the same protein, setting `merge=True` will group them into one
category and the resulting metrics will be averages of the grouped copies. Setting `merge=False` will treat each protein individually.
If you are providing a homopolymer protein, then you can use this option to keep the chains separated.

#### Using topologies

```python
# Get the lipids defined in the system:
lipids.lipid_names()
> array(['POPE', 'POPS', 'CHOL'], dtype=object)  # example output

# Get the lipid name and count:
lipids.lipid_count()
> {'POPE': 652, 'POPS': 652, 'CHOL': 652}  # example output

# Get the lipid indices:
lipids.l_indices
< array([10780, 10781, 10782, ..., 23817, 23818, 23819], dtype=int64)  # example output

# proteins is a list:
print (proteins)
> [<prolintpy.Protein containing 4 replicate(s) of Protein0 and 273 beads each>]  # example output

# Number of protein copies:
protein = proteins[0]
protein.counter
> 4  # example output
```

```python
# Describe the protein topology:
print (f'The input system provided to prolintpy contains {len(proteins)} protein(s). \
\n{"~".join([f"Protein number {i+1} contains {x.counter} copies. Each copy has {protein.n_residues} residues: {protein.first_residue}-{protein.last_residue} range and {len(protein.beads)} atoms/beads." for i, x in enumerate(proteins)])}\
'.replace('~', '\n'))

  # example output
> The input system provided to prolintpy contains 4 protein(s).
  Protein number 1 contains 1 copies. Each copy has 273 residues: 1-273 range and 630 atoms/beads.
  Protein number 2 contains 1 copies. Each copy has 273 residues: 1-273 range and 630 atoms/beads.
  Protein number 3 contains 1 copies. Each copy has 273 residues: 1-273 range and 630 atoms/beads.
  Protein number 4 contains 1 copies. Each copy has 273 residues: 1-273 range and 630 atoms/beads.
```

### Customizing definitions

Using the code above we created topologies for the whole system, with a very simple syntax. We can also customize our input.
For instance, if our simulations are at the atomistic level of detail, all we have to do is define: `resolution="atomistic"` and prolintpy will
work the same way for atomistic data. If we want to only focus on one or a subsample of lipids (instead of all lipids in the system) we can do that
by defining the `lipid_names` option, like so:

```python
# Define topologies for an atomistic system, but only considering cholesterol lipids and keeping chains separate:
resolution = "atomistic"
lipids = pl.Lipids(t.topology, resolution=resolution, lipid_names=['CHL1'])
proteins = pl.Proteins(t.topology, resolution=resolution).system_proteins(merge=False)
```

There is also an `atom_names` option that you can supply to `pl.Lipids` to futher narrow down the lipids you want to
consider by only focusing on a specific atom/bead of a particular lipid or lipid type.

It would be possible to combine the lipid and protein topologies into one description, however, the decision to have separate topologies
was done for better readability of the code and future proofing to allow any type of lipid (or ligand) to be represented.
