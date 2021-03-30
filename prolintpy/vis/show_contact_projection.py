import numpy as np
import mdtraj as md
import matplotlib as mpl
from matplotlib.pyplot import cm

# import nglview as nv

from prolintpy.utils.shift_range import shift_range

def show_contact_projection(t, bf, protein, offset=1, ngl_repr='surface', cmap='Reds', only_backbone=True):
    """Visualize lipid-protein contacts by mapping them onto the structure of the protein.

    Parameters
    ----------

    t : MDTraj.Trajectory

    bf : list
        List of contacts.

    protein : ProLint.Protein

    offset : int
        Useful when the protein residue numbering does not start from 1.
        Default is 1.

    ngl_repr: str
        One representation that will be used by nglview to display the protein. The following are supported
        options: 'point', 'line', 'rope', 'tube', 'trace', 'label', 'cartoon', 'licorice', 'ribbon',
        'backbone', 'spacefill', 'ball+stick'

    cmap : str
        One of matplotlib cmaps to be used to color contacts.

    only_backbone : bool
        Display only the backbone atoms. Only True is currently supported. Default is True.


    """

    df = protein.dataframe[0]

    # only_backbone = True,
    if only_backbone:
        if protein.resolution == "martini":
            indices = df[(df.name == "BB")].index.to_numpy()
        elif protein.resolution == "atomistic":
            indices = df[(df.name == "CA")].index.to_numpy()
    else:
        indices = df.index.to_numpy()

        ind_array = protein.get_indices()
        b = []
        for i, bfac in enumerate(bf):
            for rep in ind_array[i]:
                b.append(bfac)

        bf = b


    t_slice = t[0].atom_slice(indices)


    bf_cmap = cm.get_cmap(cmap)

    colors = [mpl.colors.to_hex(x) for x in bf_cmap(shift_range(bf))]
    cs = [[y, str(x+offset)] for x, y in enumerate(colors)]

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


