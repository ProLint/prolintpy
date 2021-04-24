import numpy as np
import pandas as pd
import mdtraj as md
from collections import defaultdict, Counter
from prolintpy.utils.martinilipids import martini_lipids, lipid_species, martini_headgroups
from prolintpy.utils.metrics import FUNCS, calculate_metric, residence_time
# from prolintpy.core.systemtopology import Lipids

class ComputeContacts(object):
    """Stores information on the contact analysis between system proteins and lipids.

    Instantiating the class only creates the object and populates a few attributes. To calculate
    and analyze contacts you have to call the compute_neighbors method which is the most time
    consuming part and depending on the input parameters can take considerable time to finish.

    Attributes
    ----------
    t : MDTraj.Trajectory

    proteins : ProLint.Proteins.system_proteins()

    lipids : ProLint.Lipids

    """

    def __init__(self, t, proteins, lipids, timestep=None):
        """Create new ProLint Contacts class"""

        """

        timestep is used to convert contacts into trajectory time units. The simplest case,
        when stride=1, is to multiply contacts with dt. When a stride > 1 is used, however, the
        exact time cannot be retrieved. ProLint can estimate that using either:
        n_frames * stride - stride, or using (n_frames * stride) - (stride - 1).

        """

        # TODO:
        # assert that proteins is a list

        self.frames = (t.n_frames - 1)
        self.timestep = t.timestep
        self.t_unit = "unit"
        self.dt = None

        self.proteins = proteins
        self.lipids = lipids

    def compute_neighbors(self, t, residues=None, cutoff=0.7, proteins=None, atom_names=[], grouped=True, mdtraj=False, leightweight=False):
        """Compute contacts between proteins and given lipids in the system.

        Parameters
        ----------
        t : MDTraj.Trajectory
        residues : list
            Explicitly specify the residues for contact calculations. Really important to speed
            up calculations when you are only interested in a small part of the protein.
            It should be 1-based, meaning the residue numbers should correspond to residues in the input topology.
        cutoff : float
            The maximum distance that two atoms/beads are considered to be in contact.
        atom_names : list
            Restrict selection of lipids only to the following atoms/beads.
        grouped : bool, default=True
            Group lipids in the system according to their headgroup type. E.g. POPC and DOPC
            would be grouped into PC lipids, etc.
        mdtraj : bool, default=False
            Store MDTraj calculated contacts.

        Returns
        -------
        prolint_contacts : dict
            Dictionary of ProLint contacts. Contacts for each residue of each protein
            are stored using the ProLint.LPContacts class.

        """

        if len(atom_names) == 0:
            haystack_indices = self.lipids.l_indices
            lipid_names = self.lipids.lipid_names()
        else:
            df = self.lipids.dataframe
            haystack_indices = list(df[df.name.isin(atom_names)].index)
            lipid_names = df[df.name.isin(atom_names)].resName.unique()

        if grouped:
            PLASMA_LIPIDS = martini_lipids(lipid_names)
        else:
            PLASMA_LIPIDS = {}
            for lip in lipid_names:
                PLASMA_LIPIDS[lip] = [lip]

        if proteins is None:
            proteins = self.proteins

        prolint_contacts = defaultdict(dict)
        for protein in proteins:
            for pc in np.arange(protein.counter):
                print ("Working on protein copy: {}".format(pc))

                if residues is None:
                    residues = protein.residues

                residue_indices = protein.get_indices(residues=residues, df=protein.dataframe[pc])
                per_residue_results = {}
                for chain in range(1, protein.chains+1):

                    for r, idxs in enumerate(residue_indices):
                        # print (residues[r])

                        atoms_count = int(idxs.size / protein.chains)
                        idx = idxs[atoms_count*(chain-1):atoms_count*chain]

                        chain_residue = residues[r] + (protein.n_residues * (chain-1))

                        md_out = md.compute_neighbors(t[1::], cutoff, idx, haystack_indices=haystack_indices)
                        per_residue_results[chain_residue] = LPContacts(t.topology, md_out, PLASMA_LIPIDS,
                                                                self.frames, self.timestep, residue=residues[r])

                prolint_contacts[protein.name][pc] = per_residue_results

        return dict(prolint_contacts)


    def compute_lipid_distances(self, t, protein, distances_dict, PLASMA_LIPIDS, lipids_found, resolution="martini", p=False, grouped=True):

        import itertools

        def get_lipid_atoms(lipid, ldf):
            atoms = ldf[ldf.resName == lipid].name.unique().tolist()
            # getting attoms in this order: P, Ox, Cx, random
            if 'P' in atoms:
                return 'P'
            elif 'O' in ''.join(atoms):
                o = [x for x in atoms if x.startswith('O')]
                o.sort()
                return o[0]
            elif 'C' in ''.join(atoms):
                c = [x for x in atoms if x.startswith('C')]
                c.sort()
                return c[0]
            else:
                return np.random.choice(atoms)

        if p:
            ldf = p.ldf
        else:
            ldf = t.topology.to_dataframe()[0]
        lipid_residue_distances = {}
        lipids = list(distances_dict.keys())

        for lipid in lipids:
            if resolution == "martini":
                headgroup = martini_headgroups(lipid, False)[lipid]
            else:
                headgroup = get_lipid_atoms(lipid, ldf)

            for residue in distances_dict[lipid]:
                residue_distances = {}
                for prot_idx in range(protein.counter):
                    indices = protein.get_indices([residue], copy=prot_idx, suppress=True)[0]

                    df_res = protein.dataframe[prot_idx].loc[indices]
                    if protein.resolution == "martini":
                        residue_atom = "BB"
                    elif protein.resolution == "atomistic":
                        residue_atom = "CA"

                    bb = df_res[df_res.name == residue_atom].index.to_list()
                    if len(bb) == 0:
                        continue

                    lipid_idxs = ldf[
                        (ldf.name == headgroup) &
                        (ldf.resName == lipid)
                    ].index.to_list()

                    d = md.compute_distances(t, list(itertools.product(bb, lipid_idxs)))

                    dist_dict_values = {}
                    dist_sum_list = []
                    dist_label_list = []
                    for dist in range(d.shape[1]):
                        sort_array = d[:, dist].copy()
                        sort_array.sort()

                        dist_sum_list.append(sum(sort_array[-int(len(sort_array)*0.1):]))
                        dist_label_list.append(lipid_idxs[dist])
                        dist_dict_values[lipid_idxs[dist]] = d[:, dist]

                    strongest_binding_label = dist_label_list[np.argmin(dist_sum_list)]
                    residue_distances[prot_idx] = dist_dict_values[strongest_binding_label]

                if lipid_residue_distances.get(lipid):
                    lipid_residue_distances[lipid][residue] = residue_distances
                else:
                    lipid_residue_distances[lipid] = {residue: residue_distances}

        tb = {}
        for k, v in PLASMA_LIPIDS.items():
            if k in lipids_found:
                tightest_lipid = {}
                for vv in v:
                    for resi in lipid_residue_distances[vv].keys():
                        sum_value = []
                        for x in lipid_residue_distances[vv][resi].values():
                            lower_best = np.sort(x)[:int(len(x) * 0.2)]

                            sum_value.append(sum(lower_best))
                        sum_value = sum(sum_value)
                        if tightest_lipid.get(resi):
                            tightest_lipid[resi].append(sum_value)
                        else:
                            tightest_lipid[resi] = [sum_value]

                for kk, sums in tightest_lipid.items():
                    lipid_bound = v[np.argmin(sums)]
                    if tb.get(k):
                        tb[k][kk] = lipid_bound
                    else:
                        tb[k] = {kk: lipid_bound}

        final_results = {}
        for k, v in tb.items():
            for resk, resv in v.items():
                if final_results.get(resv):
                    final_results[resv][resk] = lipid_residue_distances[resv][resk]
                else:
                    final_results[resv] = {resk: lipid_residue_distances[resv][resk]}

        return final_results, dict(time=t.time, protein=[protein.name])


    def compute_distances(self, t, protein, residues, lipid, headgroup, residue_atom=None, distance_co=0.7, percentile_co=0.05, n=None):
        """Compute the distance between a protein residue and lipid atom/bead.

        Parameters
        ----------
        t : MDTraj.Trajectory
            An MDTraj.Trajectory instance.

        protein : ProLint.Protein
            An ProLint.Protein instance.

        residues : int
            The residue number for distance calculation. Numbering starts at 1.

        lipid : str
            The name of the lipid for distance calcualtion.

        headgroup : str
            The name of the atom/bead of the lipid for distance calcualtion. This is usually the headgroup.

        residue_atom : str
            The name of the atom/bead of the residue. This is by default the carbon alpha or backbone bead.
            You can specify any residue atom/bead however.

        distance_co : float
            A cutoff distance that a lipid must satisfy. This cutoff is very useful to ensure not all distances are displayed.
            Only lipids that satisfy this cutoff for percentile_co time are stored.

        percentile_co : float
            The percentage of the trajectory (measured in frames) that a lipid must be within the cutoff for it to be stored. Default is 0.05 (5%) which means
            only lipids that are at least 5% of the time within the default cutoff (0.7) are stored.

        n : dict
            ProLint contact dictionary. If this has already been calculated, it can be supplied to speed up the calculation of
            distances because only lipids that are known to have interacted with the residue are checked. It is optional however.
        """

        import itertools

        whole_dist_dict = {}

        for residue in residues:
            dist_dict = {'Time' : t.time}
            column_count = {}
            min_data = []
            for prot_idx in range(protein.counter):
                indices = protein.get_indices([residue], copy=prot_idx, suppress=True)[0]

                df_res = protein.dataframe[prot_idx].loc[indices]
                if residue_atom is None:
                    if protein.resolution == "martini":
                        residue_atom = "BB"
                    elif protein.resolution == "atomistic":
                        residue_atom = "CA"

                bb = df_res[df_res.name == residue_atom].index.to_list()
                if len(bb) == 0:
                    continue

                ldf = t.topology.to_dataframe()[0]
                if n is None:
                    lipid_idxs = ldf[
                        (ldf.name == headgroup) &
                        (ldf.resName == lipid)
                    ].index.to_list()

                else:
                    lipid_idxs = []
                    lipid_resIDs = n[protein.name][prot_idx][residue].lipids[lipid]
                    if len(lipid_resIDs) == 0:
                        continue
                    for resid in lipid_resIDs:
                        lipid_idx = ldf[
                            (ldf.name == headgroup) &
                            (ldf.resSeq == resid)
                        ].index.to_list()[0]

                        lipid_idxs.append(lipid_idx)

                d = md.compute_distances(t, list(itertools.product(bb, lipid_idxs)))

                column_count[prot_idx] = 0
                for dist in range(d.shape[1]):
                    sort_array = d[:, dist].copy()
                    sort_array.sort()
                    max_allowed = sort_array[int(len(sort_array)*percentile_co)]

                    if max_allowed < distance_co:
                        dist_dict[f'{prot_idx}_{lipid_idxs[dist]}'] = d[:, dist]
                        column_count[prot_idx] += 1
                        min_data.append(max_allowed)

            if len(min_data) == 0:
                min_val = 0
            else:
                min_val = min(min_data)

            meta = dict(
                protein=protein.name,
                lipid=lipid,
                residue=residue,
                headgroup=headgroup,
                residue_atom=residue_atom,
                cutoff=0.7,
                protein_id = list(range(protein.counter)),
                column_count = column_count,
                min_val = min_val
                )

            whole_dist_dict[residue] = (pd.DataFrame(dist_dict), meta)

        return whole_dist_dict



class LPContacts(object):
    """Stores lipid-protein contact information in a readible and usable format.

    Processes the output of the MDTraj.compute_neighbors(), calculating contacts and
    outputing the result in a readable object.

    Attributes
    ----------
    topology : MDTraj.Topology

    mdtraj_contacts : MDTraj.compute_neighbors()

    PLASMA_LIPIDS : dict
        ProLint lipids grouped according to the grouped attribute
        of ProLint.ComputeContacts.compute_neighbors(). E.g. dict(PC=["POPC, "DOPC"]) if grouped=True or
        dict(POPC=["POPC"], DOPC=["DOPC"]) if grouped=False.

    frames : int
        E.g. MDTraj.Trajectory.n_frames

    timestep : int

    residue : int or str
    """
    # Requires import of lipid_species
    def __init__(self, topology, mdtraj_contacts, PLASMA_LIPIDS, frames, timestep, mdtraj=False, residue='residue'):
        """Create LPContacts class and populate attributes"""

        all_lipid_name_id = {}
        occupancy_binaries = {key: np.zeros(frames) for key in PLASMA_LIPIDS.keys()}

        for index, values in enumerate(mdtraj_contacts):

            lipid_ids = {}
            for value in values:

                bead = topology.atom(int(value))
                lipid_name = lipid_species(bead.residue.name, PLASMA_LIPIDS)

                lipid_ids[bead.residue.resSeq] = lipid_name
                occupancy_binaries[lipid_name][index] = 1

            for k, v in lipid_ids.items():
                if all_lipid_name_id.get(v):
                    all_lipid_name_id[v].append(k)
                else:
                    all_lipid_name_id[v] = [k]

        contacts, lipids = {}, {}
        for k in PLASMA_LIPIDS.keys():
            if all_lipid_name_id.get(k):
                contact_list = list((Counter(all_lipid_name_id[k]).values()))
                contacts[k] = [x*timestep for x in contact_list]
                lipids[k] = np.unique(all_lipid_name_id[k])
            else:
                contacts[k] = [0]
                lipids[k] = []

        self.contacts = contacts
        self.occupancy = occupancy_binaries
        self.residue = residue
        self.lipids = lipids

        if mdtraj:
            self.raw_data = mdtraj_contacts


    def __str__(self):
        return "<prolintpy.LPContacts for residue {}>".format(self.residue)


    def __repr__(self):
        return "<prolintpy.LPContacts for residue {}>".format(self.residue)



def retrieve_contacts(n, r, protein=None, contacts='contacts'):
    """Get contact information for a given residue. A similar function is
    retrieve_contacts_flat.

        Parameters
        ----------
        n : ProLint.ComputeContacts
            The dictionary output of ProLint contact calculations.

        r : int
            The residue number (not index)
        protein: str
            If more than one protein in the system, you can specify the protein using
            its name.
        contacts : str
            Type of contact to calculate. 'contacts' or 'occupancy' are allowed.

        Returns
        -------
        contacts : dict
            dictionary of key: value pairs representing lipid: contact * nr_replicates

    """

    if (len(n.keys()) != 1) & (protein is None):
        raise ValueError("""prolintpy.Contacts provided contains information for more than 1 protein. \
Indicate for which protein to calculate contacts using the 'protein' option.""")
    elif protein is None:
        protein = list(n.keys())[0]

    contact_data = defaultdict(list)

    # Get the first key of the input dictionary.
    key1 = list(n.keys())[0]
    key2 = 0
    key3 = list(n[key1][0].keys())[key2]

    is_class = n[key1][key2][key3]

    if isinstance(is_class, LPContacts):
        for pc in n[protein].keys():
            if contacts == 'contacts':
                contact = n[protein][pc][r].contacts
            elif contacts == 'occupancy':
                contact = n[protein][pc][r].occupancy
            else:
                raise ValueError("Only contacts and occupancy are currently implemented.")

            for lipid, values in contact.items():

                contact_data[lipid].append(np.array(values))
    else:
        for pc in n[protein].keys():
            if contacts == 'contacts':
                contact = n[protein][pc][r]
            else:
                raise ValueError("When leightweight is used, only contacts are calcualted (no occupancies, etc.).")

            for lipid, values in contact.items():
                contact_data[lipid].append(np.array(values))

    return contact_data


def retrieve_contacts_flat(n, r, lipid=None, protein=None, contacts='contacts'):
    """Get contact information for a given residue. A similar function is
    retrieve_contacts.

        Parameters
        ----------
        n : ProLint.ComputeContacts
            The dictionary output of ProLint contact calculations.

        r : int
            The residue number (not index)
        lipid : list
            A list of lipids for which to calculate contacts.
        protein: str
            If more than one protein in the system, you can specify the protein using
            its name.
        contacts : str
            Type of contact to calculate. 'contacts' or 'occupancy' are allowed.

        Returns
        -------
        contacts : array
            1D array of all contacts for all protein copies in the system.

    """

    if (len(n.keys()) != 1) & (protein is None):
        raise ValueError("""prolintpy.Contacts provided contains information for more than 1 protein. \
Indicate for which protein to calculate contacts using the 'protein' option.""")
    elif protein is None:
        protein = list(n.keys())[0]

    def is_lipid_defined(c, lipid):
        if (len(c.keys()) != 1) & (lipid is None):
                raise ValueError("""prolintpy.Contacts provided contains information for more than 1 lipid. \n
Indicate for which lipid to calculate contacts using the 'lipid' option.""")
        elif lipid is None:
            lipid = list(c.keys())[0]

        return lipid



    # Get the first key of the input dictionary.
    key1 = list(n.keys())[0]
    key2 = 0
    key3 = list(n[key1][0].keys())[key2]

    is_class = n[key1][key2][key3]

    contact_data = defaultdict(list)
    if isinstance(is_class, LPContacts):
        test_contact = n[protein][0][r].contacts
        lipid = is_lipid_defined(test_contact, lipid)

        for pc in n[protein].keys():
            if contacts == 'contacts':
                contact = n[protein][pc][r].contacts
            elif contacts == 'occupancy':
                contact = n[protein][pc][r].occupancy
            else:
                raise ValueError("Only contacts and occupancy are currently implemented.")

            for l, v in contact.items():
                if l != lipid:
                    continue
                contact_data[l].append(np.array(v))

    else:
        test_contact = n[protein][0][r]

        lipid = is_lipid_defined(test_contact, lipid)

        for pc in n[protein].keys():
            if contacts == 'contacts':
                contact = n[protein][pc][r]
            else:
                raise ValueError("When leightweight is used, only contacts are calcualted (no occupancies, etc.).")

            for l, v in contact.items():
                if l != lipid:
                    continue
                contact_data[l].append(np.array(v))

    return np.concatenate(contact_data[lipid])




def contacts_dataframe(n, p, t, radius, resolution="martini", output_erros=False, co=0, custom_metrics=None, time_unit='us', rt=False, range_type='geo', norm=True):
    """Convert contact information into a DataFrame format. The function calculates
    several contact metrics for each residue of each unique protein in the system.
    Contact information for multiple copies of the same protein are averaged. The following
    are the contacts supported:

    ``mean`` (Mean_Duration) : The average duration of all contacts.

    ``max`` (Longest_Duration) : The longest duration (averaged if more than 1 protein).

    ``sum`` (Sum_of_all_Contacts) : The total sum of all contacts.

    ``lnr`` (Lipid_Number) : The average number of lipids in contact (the total number of contacts
    normalized with respect to the number of frames).

    ``nlnr`` (Normalized_Lipid_Number) : Lipid_Number normalized with respect to the
    number of different lipids (e.g. number of different cholesterols).

    ``occ`` (Occupancy) : For each frame, we give a value of 0 if no lipid of interest is in
    contact with the residue and 1 otherwise. Occupancy is then: sum_of_all_1s / nr_of_frames.

    ``rt``  (Residence_Time) : See here: https://pubs.acs.org/doi/abs/10.1021/ja310577u
    and here: https://www.pnas.org/content/117/14/7803

    Residence time is quite time consuming and by default is not calculated. This can be altered
    using the rt option.

    Parameters
    ----------
    n : ProLint.ComputeContacts
        The dictionary output of ProLint contact calculations.

    p : ProLint.Proteins.system_proteins()
        list containing protein information for the system.

    t : MDTraj.Trajectory

    radius : float
        the cutoff radius used when calculating contacts.

    resolution: string
        the resolution of the data. Default martini.

    output_erros: bool
        Calculate errors and output them in the dataframe? Default is False.
        Currently, errors are only meaningfull if there are multiple proteins in the system.

    co : int
        Discard contacts shorter than co value. Useful when you want to discard contacts
        with a short duration.

    custom_metrics : list
        Provide a custom list of parameters (from the default/supported ones). E.g.
        custom_metrics=['occ'] to calculate only the occupancy.

    time_unit : str
        The timeunit of contact measurements. Either 'us' for microsecond or 'ns' for nanosecond.
        Only used if norm=True.

    rt : bool, default=False
        Calculate the residence time.

    range_type : str
        Residence time can be sped up by supplying a default range. You can specify either a geometric,
        linear, or mixed range ('geo', 'mixed', 'linear'). Fastest and default is 'geo'.

    norm : bool, default=True
        Normalize contacts with respect to time_unit.

    Returns
    -------
    dataframe : pandas.DataFrame
        DataFrame containing contact information for the default metrics.

    """

    RESULTS = defaultdict(list)
    metrics = {
        'mean':'Mean_Duration',
        'max' :'Longest_Duration',
        'sum' :'Sum_of_all_Contacts',
        'lnr' :'Lipid_Number',
        'nlnr':'Normalized_Lipid_Number',
        'occ' :'Occupancy',
    }
    if rt:
        metrics['rt'] = 'Residence_Time'

    if custom_metrics is not None:
        metrics = {key : metrics[key] for key in custom_metrics}

    if resolution == "martini":
        ref_atom = 'BB'
    elif resolution == "atomistic":
        ref_atom = 'CA'

    proteins = list(n.keys())
    for prot_idx, protein in enumerate(proteins):
        residues = list(n[protein][0].keys())

        # Get residue names based on the indices of unique resSeq occurrances.
        # This should work for any ff supported by mdtraj not just cg.
        df = p[prot_idx].dataframe[0]
        resnames = df[df.name == ref_atom].resName.to_numpy()

        for res_idx, residue in enumerate(residues):
            for k, v in metrics.items():
                if k == 'occ':
                    con = retrieve_contacts(n, residue, protein=protein, contacts='occupancy')
                else:
                    con = retrieve_contacts(n, residue, protein=protein)

                if k != 'rt':
                    if (k in ['mean', 'max', 'sum']) & (norm):
                        result = calculate_metric(con, FUNCS[k], co, t, norm=norm, unit=time_unit)
                    else:
                        result = calculate_metric(con, FUNCS[k], co, t)
                else:
                    result = {}
                    for i, j in con.items():
                        sdf = t.topology.to_dataframe()[0]
                        ln = sdf[sdf.resName == i].resSeq.unique().size
                        time = t.time[-1] - t.time[0]
                        result[i] = residence_time(np.concatenate(j), ln, time, range_type=range_type)

                for lipid, values in result.items():
                    RESULTS[v].append(values[0])
                    if output_erros:
                        RESULTS[v + '_Error'].append(values[1])
                    if k == list(metrics.keys())[-1]:
                        RESULTS["Protein"].append(protein)
                        RESULTS["Lipids"].append(lipid)
                        RESULTS["Radius"].append(float(radius))
                        RESULTS["ResID"].append(residue)
                        RESULTS['ResName'].append(resnames[res_idx])

    return pd.DataFrame(RESULTS).sort_values(['Protein', 'Lipids', 'ResID']).reset_index(drop=True)



def retrieve_distances(protein_dataframe, group_lipids, resolution, lipids, top_nr=30):

    if group_lipids and resolution == "martini":
        SYSTEM_LIPIDS = martini_lipids(lipids.lipid_names())
    else:
        SYSTEM_LIPIDS = {}
        for lip in lipids.lipid_names():
            SYSTEM_LIPIDS[lip] = [lip]


    distances_dict = {}
    for protein in protein_dataframe.Protein.unique():
        df = protein_dataframe[protein_dataframe.Protein == protein]
        idx = df.Longest_Duration.sort_values(ascending=False)[:int(top_nr)].index
        df = df[df.index.isin(idx)]
        lipids_found = df.Lipids.unique()
        for lipid in lipids_found:
            if group_lipids:
                for group in SYSTEM_LIPIDS[lipid]:
                    distances_dict[group] = df[df.Lipids == lipid].ResID.to_list()
            else:
                distances_dict[lipid] = df[df.Lipids == lipid].ResID.to_list()

    return distances_dict, SYSTEM_LIPIDS, lipids_found
