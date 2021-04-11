import re
import numpy as np
import pandas as pd

class SystemTopology(object):
    """Base class for Lipids and Proteins classes to inherit from.

    It reads an MDTraj topology and extracts useful information. Its main
    purpose is intended to be in distinguishing the resolution of the input
    system. Other ProLint classes and functions should ideally be
    resolution-agnostic, though currently that is not entirely the case.

    ProLint, in its current version, only supports the Martini model.


    Attributes
    ----------
    topology : MDTraj.Topology

    resolution : Either "martini" or "atomistic".

    """

    def __init__(self, topology, resolution="martini"):
        """Create a new ProLint topology."""

        if resolution == "martini":
            lipid_indices = topology.select("all and not (name BB or (name =~ 'SC[1-9]'))")
            proteins_indices = topology.select("name BB or (name BB or (name =~ 'SC[1-9]'))")
        elif resolution == "atomistic":
            lipid_indices = topology.select("all and not protein")
            proteins_indices = topology.select("protein")
        else:
            raise NotImplementedError('Provided resolution is not supported. Only "martini" and "atomistic" are supported.')

        self.resolution = resolution
        self.topology = topology
        self.dataframe = topology.to_dataframe()[0]
        self.l_indices = lipid_indices
        self.p_indices = proteins_indices
        self.pdf = self.dataframe.iloc[proteins_indices]
        self.ldf = self.dataframe.iloc[lipid_indices]


    def __str__(self):
        return "Base Class for the topology of the MDTraj system."

    def __repr__(self):
        return "Base Class for the topology of the MDTraj system."


class Lipids(SystemTopology):
    """
    Lipids stores information about the lipid composition of the system.
    Plasma membrane lipids (as described here: https://pubs.acs.org/doi/abs/10.1021/acscentsci.8b00143)
    work best, but all lipids in the Martin lipidome should work (see here:
    http://www.cgmartini.nl/index.php/force-field-parameters/lipids)

    g_surf only works with Plasma lipids.

    Attributes
    ----------
    topology : MDTraj.Topology
        MDTraj.Topology instance

    lipid_names : list
        Only select the following lipids. Useful when you want to analyze contacts
        with only one or a few lipids instead of all system lipids.

    atom_names : list
        Only select the following atoms/beads. Useful when you want to further
        restrict analysis to specific bead/atom types instead of whole lipid residue.

    """

    def __init__(self, topology, resolution="martini", lipid_names=None, atom_names=None):
        """Create a new ProLint Lipids class."""

        super().__init__(topology, resolution=resolution)

        lfunc = lambda x, y: True if y is None else x in y
        user_indices = [a.index for a in topology.atoms
                        if lfunc(a.residue.name, lipid_names)
                        and lfunc(a.name, atom_names)
                        ]

        if len(user_indices) == 0:
            print ("No lipids have been selected!.")

        self.l_indices = np.intersect1d(self.ldf.index.to_numpy(), user_indices)
        self.dataframe = self.dataframe.iloc[self.l_indices]


    def lipid_names(self):
        """Get the names of all lipids that will be analyzed.

        Returns
        -------
        names : array of lipid names
        """
        return self.dataframe.resName.unique()

    def lipid_count(self):
        """Get the name and count of each lipid that will be analyzed.

        Returns
        -------
        lipic_count : dictionary
            key:value corresponds to lipid_name:count.
        """
        lc = {}
        lipids = self.lipid_names()
        for lipid in lipids:
            lc[lipid] = self.dataframe[self.dataframe.resName == lipid].resSeq.unique().size
        return lc


    def __str__(self):
        return "<prolintpy.Lipids containing {} atom(s)>".format(self.l_indices.size)

    def __repr__(self):
        return "<prolintpy.Lipids containing {} atom(s)>".format(self.l_indices.size)


class Protein(object):
    """
    Store various information for each protein in the system.

    Attributes
    ----------
    name : str
        The name of the protein.

    """


    def __init__(self, name):
        """Create a new ProLint Protein class."""
        self.name = name
        self.counter = 1
        self.dataframe = []
        self.beads = 0
        self.n_residues = 0
        self.residues = None
        self.resnames = None
        self.chains = 1
        self.resolution = None

    def get_indices(self, residues=None, df=None, copy=None, resolution='martini', suppress=False):
        """Return the indices for each protein residue. If there are more than one
        protein replicates/copies in the system the method will choose the residue indices of
        the first occuring protein. Alternatively you can provide the dataframe of the protein
        copy yourelf using the 'dataframe' attribute of this class.

        Parameters
        ----------
        residues : list
            Explicitly specify the residues. List should be 1-based, meaning the residue numbers should
            correspond to residues in the input topology.
        df : DataFrame
            The protein dataframe. It can be a custom definition, but usually you would want
            one of the dataframes contined in the self.dataframe attribute.
        copy : int
            Specify the protein copy/replicate.
        resolution : str
            The resolution of the system. Either "martini" or "atomistic".
        suprress : bool, default=True
            Suppress warning messages.

        Returns
        -------
        indices : list
            list of numpy arrays.
        """
        if df is None:
            if not suppress:
                if len(self.dataframe) > 1 and copy is None:
                    print (re.sub(' +', ' ', """Object has more than one dataframe. Using the first one.
                            Alternatively, you can specify the dataframe of the protein yourself using the
                            dataframe attribute of this class."""))
                else:
                    print ("Using the available dataframe")
            df=self.dataframe[0]

        if copy is not None:
            df=self.dataframe[copy]

        if residues is None:
            if resolution == "martini":
                residues = df[df.name == "BB"].resSeq.to_list()
            elif resolution == "atomistic":
                residues = df[df.name == "CA"].resSeq.to_list()
            else:
                raise ValueError("Only martini and atomistic are currently supported.")

            residues = residues[:int(len(residues) / self.chains)]
        else:
            assert all(x in self.residues for x in residues), "Supplied residues should correspond to input topology residues."

        # if self.chains > 1:
        #     residues = residues[:int(len(residues) / self.chains)]

        indices = [df[df.resSeq == x].index.to_numpy() for x in residues]

        return indices

    def __str__(self):
        return "<prolintpy.Protein containing {} replicate(s) of {} and {} beads each>".format(self.counter, self.name, self.n_residues)

    def __repr__(self):
        return "<prolintpy.Protein containing {} replicate(s) of {} and {} beads each>".format(self.counter, self.name, self.n_residues)


class Proteins(SystemTopology):
    """
    Proteins stores information about the protein composition of the system.

    Gromacs coordinate files do not contain information on protein names and their
    count. This class tries to calculate and extract this information from the topology.

    The class has a 'system_proteins' method that can be used to retrieve information
    about each protein specifically.

    Attributes
    ----------
    topology : MDTraj.Topology
        An MDTraj.Topology instance.

    """

    def __init__(self, topology, resolution="martini"):
        """Create a new ProLint Proteins class."""

        super().__init__(topology, resolution=resolution)

        # Get start and end indices of proteins in the system.
        # The assumption here is that proteins are ordered and the start residue of the next
        # protein is always smaller than the last residue of the previous protein.
        resseq = self.pdf.resSeq.to_list()
        p0 = resseq[0]
        # system_proteins contains the start and end indices of all proteins.
        fi_li = []
        fi = 0
        for li, p in enumerate(resseq):
            if p < p0:
                fi_li.append((fi, li-1))
                fi = li
            p0 = p
        fi_li.append((fi, li))

        self.fi_li = fi_li

    def system_proteins(self, merge=True):
        """Stores information for each protein in the system using the Protein class.

        Returns
        -------
        proteins : list
            list of Protein classes. One for each different protein in the system.

        """
        if self.resolution == "martini":
            ref_atom = 'BB'
        elif self.resolution == "atomistic":
            ref_atom = 'CA'

        c = 0
        proteins = []
        # Two proteins are the same if residue number and all beads are equal between them.
        for values in self.fi_li:

            # first and last index
            fi = values[0]
            li = values[1]

            # Get the data for the current protein.
            current_protein_df = self.pdf[(self.pdf.index >= fi) & (self.pdf.index <= li)]
            curr_len = len(current_protein_df[current_protein_df.name == ref_atom])

            # curr_len = len(current_protein_df)
            curr_names = current_protein_df.name.to_list()

            new_protein = True
            if merge:
                for pc in proteins:
                    if curr_names == pc.beads:
                        pc.counter += 1
                        pc.dataframe.append(current_protein_df.copy())
                        new_protein = False

            if new_protein:
                protein = Protein(f'Protein{c}')
                protein.dataframe.append(current_protein_df.copy())
                protein.beads = curr_names
                protein.n_residues = curr_len
                protein.residues = current_protein_df.resSeq.unique()

                resnames = current_protein_df[current_protein_df.name == ref_atom].resName.to_numpy()
                protein.resnames = resnames

                protein.first_residue = current_protein_df.resSeq.to_list()[0]
                protein.last_residue = current_protein_df.resSeq.to_list()[-1]

                protein.resolution = self.resolution

                proteins.append(protein)

            c += 1

        return proteins

    def merge_chains(self, proteins):
        """If the input structure file contained multiple and similar chains, ProLint
        will count them as separate proteins by default. This function allows you to
        change that. You call it on the list of proteins that you want to change and it will
        merge those proteins into multiple chains.

        Parameters
        ----------
        proteins : ProLint.Proteins
            Proteins object.

        """
        merged = []
        for protein in proteins:
            protein.beads = protein.beads * protein.counter
            protein.dataframe = [pd.concat(protein.dataframe)]
            # protein.n_residues = protein.n_residues #* protein.counter => n_residues remains the same
            protein.chains = protein.counter
            protein.counter = 1
            merged.append(protein)
            print ("Merged copies as chains in {}".format(protein.name))

        return merged




    def __str__(self):
        return "class storing information on system proteins"

    def __repr__(self):
        return "class storing information on system proteins"
