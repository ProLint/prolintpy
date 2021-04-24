import os
import json
import networkx as nx
from networkx.readwrite import json_graph

from prolintpy.utils.shift_range import shift_range

def network(l, df, metric, grouped, radius=None, return_json=False):
    """Calculate network information.

    ** The code has to be updated and improved for stability and readability **

    Parameters
    ----------

    l : ProLint.Lipids
        System lipids class.

    df : pd.DataFrame
        pandas dataframe of contacts

    metric : str
        The metric to for network analysis.

    grouped : bool
        Are lipid grouped into headgroup categories.

    radius : float
        Cutoff radius used in measuring contacts. Useful when there's multiple
        cutoff radii.

    return_json : bool
        Return the results instead of saving them to a file.

    """

    RESIDUES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    if grouped:
        from prolintpy.utils.martinilipids import martini_lipids
        PLASMA_LIPIDS = martini_lipids(l.lipid_names())
    else:
        PLASMA_LIPIDS = {}
        for lip in l.lipid_names():
            PLASMA_LIPIDS[lip] = [lip]

    if radius is None:
        if df.Radius.unique().size != 1:
            raise ValueError("Ambiguous cutoff. Specify radius parameter.")
        r = df.Radius.unique()[0]
    else:
        r = radius

    ALL_LIPIDS = []
    for lip in l.lipid_names():
        lipid_df = l.ldf
        lipid_top_bead = lipid_df[(lipid_df.resName == lip)].name.unique()[0]
        c = lipid_df[(lipid_df.resName == lip) & (lipid_df.name == lipid_top_bead)].count().serial
        ALL_LIPIDS.append((lip, c))


    proteins = df.Protein.unique()
    LIPIDS = df.Lipids.unique()

    def get_data_network(r, proteins, df, metric):
        """
        Getting the data.
        """
        data = {}
        for gpcr in proteins:
            gpcr_df = df[(df.Radius == float(r)) & (df.Protein == gpcr)]
            lipids = gpcr_df.Lipids.unique()
            lipid_number = len(lipids)

            lipid_dict = {}

            for lipid in lipids:
                lipid_df = gpcr_df[gpcr_df.Lipids == lipid]
                residue_dict = {}

                for residue in RESIDUES:
                    residue_df = lipid_df[lipid_df.ResName == residue]
                    value = residue_df[residue_df[metric] > 0][metric].mean()
                    if value > 0:
                        residue_dict[residue] = value

                lipid_dict[lipid] = residue_dict

            data[gpcr] = lipid_dict


        data_dict = {}

        for gpcr in proteins:
            gpcr_df = df[(df.Radius == float(r)) & (df.Protein == gpcr)]
            lipids = gpcr_df.Lipids.unique()
            lipid_number = len(lipids)
            lipid_dict = {}

            for lipid in lipids:
                lipid_df = gpcr_df[gpcr_df.Lipids == lipid]

                res_names = lipid_df.ResName.to_list()
                res_numbers = lipid_df.ResID.to_list()
                parameters = lipid_df[metric].to_list()

                residue_dict = {}

                for res_name, res_number, parameter in zip(res_names, res_numbers, parameters):
                    name_number = res_name + "-" + str(res_number)
                    residue_dict[str(name_number)] = parameter

                lipid_dict[lipid] = residue_dict

            data_dict[gpcr] = lipid_dict

        return data, data_dict, lipid_number


    def protein_lipid_network(prot, data, data_dict, lipid_number, r, df, ALL_LIPIDS, PLASMA_LIPIDS):
        """
        Main function to calculate the protein network.
        """

        protein_df = df[(df.Radius == float(r)) & (df.Protein == prot)]
        lipids = protein_df.Lipids.unique()

        residue_id_dict = data[prot]
        protein_dict = data_dict[prot]

        residues = protein_df[protein_df.Lipids == lipids[0]].ResName
        residues_number = len(residues)
        residue_ratios = {}
        for amino_acid in RESIDUES:
            residue_number = len(residues[residues == amino_acid])
            residue_ratios[amino_acid] = residue_number/residues_number

        # this is the same as lipid_names
        system_lipids = [x[0] for x in dict(ALL_LIPIDS).items()]

        # ALL_LIPIDS may contain the same lipid twice (for each lipid).
        # When calculating the content ratio of lipids we have to consider this.
        total_lipids = 0
        lip_size = {}
        for lipid_name in system_lipids:
            leaflet_lipids = [x for x in ALL_LIPIDS if x[0] == lipid_name]
            lipid_sum = 0
            leaflet_lipid = None
            for leaflet_lipid in leaflet_lipids:
                lipid_sum += leaflet_lipid[1]
            total_lipids += lipid_sum
            if leaflet_lipid is not None:
                lip_size[leaflet_lipid[0]] = lipid_sum


        lipid_sum = 0
        lip_size_dict = {}
        for lipid_group, leaflet_lipid in PLASMA_LIPIDS.items():
            for lip in leaflet_lipid:
                if lip_size.get(lip):
                    lipid_sum += lip_size[lip]
            try:
                lip_size_dict[lipid_group] = lipid_sum/total_lipids
            except ZeroDivisionError:
                lip_size_dict[lipid_group] = 0
            lipid_sum = 0

        # We just need one lipid to slice the data frame to get the
        # residue name and number.
        residue_number_name = []
        for index in range(len(residues)):
            name = protein_df[protein_df.Lipids == lipids[0]].ResName.to_list()[index]
            number = protein_df[protein_df.Lipids == lipids[0]].ResID.to_list()[index]
            residue_number_name.append((name, number))

        # Node names and IDs
        node_number = list(lipids) + list(RESIDUES) + residue_number_name
        core_nodes = len(lipids) + len(RESIDUES)
        node_dict = {}
        count1, count2 = 0, 0
        for node in node_number:

            if count1 >= core_nodes:
                k = residue_number_name[count2]
                res_id = k[0] + "-" + str(k[1])
                node_dict[res_id] = count1
                count1 += 1
                count2 += 1
                continue

            node_dict[node] = count1
            count1 += 1


        exclusions = []
        for residue in RESIDUES:
            for name_number in residue_number_name:
                if residue == name_number[0]:

                    counter = 0
                    for lipid in lipids:
                        res_id = name_number[0] + "-" + str(name_number[1])
                        if protein_dict[lipid][res_id] > 0.05:
                            break
                        else:
                            counter += 1

                    if counter == lipid_number:
                        exclusions.append(res_id)


        node_keys, node_vals = [], []
        node_id = 0
        for node_key in node_dict.items():
            if node_key[0] not in exclusions:
                node_keys.append(node_key[0])
                node_vals.append(node_id)
                node_id += 1

        node_dict = dict(zip(node_keys, node_vals))

        # LINKING NODES
        edge_data = []
        for lipid in lipids:
            for residue in RESIDUES:
                try:
                    edge_data.append((node_dict[lipid],
                                    node_dict[residue],
                                    residue_id_dict[lipid][residue]))
                except KeyError:
                    edge_data.append((node_dict[lipid], node_dict[residue], 0))


        for residue in RESIDUES:
            for name_number in residue_number_name:
                if residue == name_number[0]:
                    for lipid in lipids:
                        res_id = name_number[0] + "-" + str(name_number[1])
                        if protein_dict[lipid][res_id] > 0.05:
                            edge_data.append((node_dict[residue], node_dict[res_id], 0.1))
                            break

        for lipid in lipids:
            for name_number in residue_number_name:
                res_id = name_number[0] + "-" + str(name_number[1])
                if protein_dict[lipid][res_id] > 0.05:
                    edge_data.append((node_dict[lipid], node_dict[res_id], protein_dict[lipid][res_id]))



        edges_tmp = []
        for edge in edge_data:

            if (edge[0] >= lipid_number) & (edge[0] < lipid_number) & (edge[1] >= lipid_number):

                for key, value in node_dict.items():
                    if value == edge[1]:
                        node_id = key

                if node_id not in exclusions:
                    edges_tmp.append(edge)

            else:
                edges_tmp.append(edge)


        # print (edges_tmp)
        weights = [x[2] for x in edges_tmp]
        weights = shift_range(weights, new_min=0.05, new_max=2)
        edges_tmp = [(x[0], x[1], weights[i]) for i, x in enumerate(edges_tmp)]
        edges = []
        for index in enumerate(edges_tmp):
            edges.append((edges_tmp[index[0]][0],
                        edges_tmp[index[0]][1],
                        {"weight": edges_tmp[index[0]][2]}))


        # Set new ranges.
        r_r = [residue_ratios[x] for x in RESIDUES]
        # node_values_1 = shift_range(min(r_r), max(r_r), r_r)
        node_values_1 = shift_range(r_r, new_min=5, new_max=15)

        l_r = [lip_size_dict[x] for x in lipids]
        # node_values_2 = shift_range(min(l_r), max(l_r), l_r)
        node_values_2 = shift_range(l_r, new_min=5, new_max=15)

        node_values = node_values_2 + node_values_1
        extend = len(node_vals) - len(node_values)
        default_size = 15 - 5
        node_values.extend([default_size] * extend)

        return (node_vals, node_values, node_keys, edges)


    data, data_dict, lipid_number = get_data_network(r, proteins, df, metric)
    collect_results = []
    for gpcr in proteins:

        results = protein_lipid_network(gpcr, data, data_dict, lipid_number, r, df, PLASMA_LIPIDS, PLASMA_LIPIDS)
        node_vals, node_values, node_keys, edges = results

        graph_nx = nx.Graph(lip_end=lipid_number-1, res_start=20+lipid_number-1)
        for i in range(len(node_vals)):
            graph_nx.add_node(i, values=node_values[i], color="red", name=node_keys[i])

        graph_nx.add_edges_from(edges)
        graph_json = json_graph.node_link_data(graph_nx)

        collect_results.append((graph_json, gpcr))

    if return_json:
        return collect_results

    for result in collect_results:
        graph_json, gpcr = result
        file_name = os.path.join(gpcr + '_' + metric + "_" + str(r) + '_network.json')
        json.dump(graph_json, open(file_name, 'w'))

    return None
