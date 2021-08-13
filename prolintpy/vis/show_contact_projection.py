import numpy as np
import mdtraj as md
import matplotlib as mpl
from matplotlib.pyplot import cm

import nglview as nv

from prolintpy.utils.shift_range import shift_range

def show_contact_projection(t, bf, protein=None, residue_list=None, ngl_repr='surface', cmap='Reds'):
    """Visualize lipid-protein contacts by mapping them onto the structure of the protein.

    Parameters
    ----------

    t : MDTraj.Trajectory

    bf : list
        List of contacts.

    protein : ProLint.Protein
    
    residue_list: list
        List of residues. Required when working with a subset of residues. 

    ngl_repr: str
        One representation that will be used by nglview to display the protein. The following are supported
        options: 'point', 'line', 'rope', 'tube', 'trace', 'label', 'cartoon', 'licorice', 'ribbon',
        'backbone', 'spacefill', 'ball+stick'

    cmap : str
        One of matplotlib cmaps to be used to color contacts.

    """

    if protein:
        first_frame = t[0]
        df = protein.dataframe[0]

        if protein.resolution == "martini":
            indices = df[(df.name == "BB")].index.to_numpy()
        elif protein.resolution == "atomistic":
            indices = df[(df.name == "CA")].index.to_numpy()

    else:
        df = t.topology.to_dataframe()[0]
        indices = df.index.to_numpy()

    t_slice = t[0].atom_slice(indices)
    resseq = df.resSeq.to_list()
    
    
    atomic_bfactors = []
    if residue_list is not None:
        init, bfactors = 1, []
        for idx, residue in enumerate(residue_list):
            while init < residue:
                bfactors.append(0)
                init += 1
            bfactors.append(bf[idx])
            init += 1
        dif = indices.size - len(bfactors)
        for remaining_residues in range(dif):
            bfactors.append(0)
            
        atomic_bfactors = []
        for atom in resseq:
                atomic_bfactors.append(bfactors[atom-1])
    else:
        if len(df.resSeq.unique()) != len(bf):
            raise TypeError ('When projecting only a subset of residues provide a list of tuples: [(residue_id, value), ...]')
        for atom in resseq:
            atomic_bfactors.append(bf[atom-1])

    bf_cmap = cm.get_cmap(cmap)

    colors = [mpl.colors.to_hex(x) for x in bf_cmap(shift_range(atomic_bfactors))]
    cs = [[y, str(resseq[x])] for x, y in enumerate(colors)]

    scheme = nv.color._ColorScheme(cs, 'bf')

    reprs = [
        'point', 'line', 'rope', 'tube', 'trace', 'label',
        'cartoon', 'licorice', 'ribbon',
        'backbone', 'spacefill', 'ball+stick'
        ]


    params = {
        'surface' : [{'type': 'surface',
                      'params': {'surfaceType': "av",
                                 'probeRadius': 2.1,
                                 'color':scheme}}],
     }
    for rep in reprs:
        params[rep] = [{'type': rep, 'params': {'color':scheme}}]

    view = nv.NGLWidget()

    view.add_component(t_slice)
    view.set_representations(params[ngl_repr])

    return view

