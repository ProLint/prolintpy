import os
import pickle
import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm

# bokeh aplication
from bokeh.document import Document

from bokeh.plotting import figure
from bokeh.layouts import layout, column, row
from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar, BasicTickFormatter
from bokeh.models.widgets import Slider, Select, RangeSlider

from bokeh.io import show, output_notebook

def flat_density(t, lipid, frames):

    if len(lipid) == 1:
        selection_string = "resname {}".format(lipid[0])
    else:

        s_st_list = ["{} or resname ".format(x) for x in lipid[:-1]]

        selection_string = "resname "
        for s in s_st_list:
            selection_string += s
        selection_string += lipid[-1]

    lipid_ndx = t.topology.select(selection_string)
    lipid_xyz = t.xyz[:, lipid_ndx, :]

    def slice_array(arr, slice_by):
        mask = np.ones_like(arr, dtype=bool)
        mask[::int(slice_by)] = False
        xshape = int(frames - mask[::int(slice_by)].shape[0])
        arr = arr[mask].reshape(xshape*arr.shape[1], 3)

        return arr

    lipid_xyz = lipid_xyz.reshape(lipid_xyz.shape[0]*lipid_xyz.shape[1], 3)

    array_memory_size = len(pickle.dumps(lipid_xyz)) / 1024
    while array_memory_size > 5000:
        lipid_xyz = slice_array(lipid_xyz, 10)

        array_memory_size = len(pickle.dumps(lipid_xyz)) / 1024

    return lipid_xyz


def show_density(t, lipids, resolution, filenames=None):
    """Visualize the preferential localization of lipids using 2D density maps.

    Parameters
    ----------

    t : MDTraj.Trajectory

    lipids: list
        List of lipids to visualize their density. Currently, only single lipids are supported. No grouping.

    resolution: string
        Either "martini" or "atomistic"

    filename: list, Optional
        A list of strings pointing to the location of npy files that are saved by ProLint.

    """

    if filenames is not None:
        lipids = []
        for filename in filenames:
            if filename.split('/')[-1].startswith('prot'):
                prot_xyz = np.load(filename)
            else:
                lipids.append(filename)
        l_xyz = np.load(lipids[0])
    else:
        from prolintpy.core.systemtopology import Proteins
        p = Proteins(t.topology, resolution=resolution)
        prot_xyz = t.xyz[0, p.p_indices, :]
        l_xyz = flat_density(t, lipids, t.n_frames)

    output_notebook()

    def densityApp(doc):
        """
        Handler function for the density application.
        """
        nonlocal t, lipids, resolution, prot_xyz, l_xyz

        all_mpl_cmaps = [
                    'viridis', 'plasma', 'inferno', 'cividis',
                    'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                    'RdYlBu', 'RdYlGn', 'Spectral','Spectral_r', 'coolwarm', 'coolwarm_r', 'seismic',
                    'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                    'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
                    'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
        ]
        lipids = [x for x in lipids if not x.startswith('prot')]

        x_min, x_max = l_xyz[:, 0].min(), l_xyz[:, 0].max()
        y_min, y_max = l_xyz[:, 1].min(), l_xyz[:, 1].max()
        z_min, z_max = l_xyz[:, 2].min(), l_xyz[:, 2].max()

        # widgets
        color_map = Select(title="Colormap", value="viridis",
                        options=all_mpl_cmaps, width=120)

        lipid = Select(title="Lipids", value=lipids[0],
                    options=lipids, width=100)

        number = Slider(title="Number of Bins", value=150, start=80, end=500, step=10, width=200)

        protein = Select(title="Show Protein", value="No",
                        options=["Yes", "No"], width=90)

        zrange = RangeSlider(start=z_min, end=z_max,
                            value=(z_min, z_max), step=0.2, title="Z-Axis Range",
                            callback_policy='mouseup', width=352)

        cbar_range = Slider(value=256, start=0, end=256, step=1, height=600,
                            orientation='vertical', callback_policy='mouseup', margin=(15, 0),
                            tooltips=False, show_value=False)

        # colorschemes
        colormap =cm.get_cmap(color_map.value)
        bokehpalette = [mpl.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]

        def lipid_density(array, bins=160):
            """Given a 2D array, containing the x, y values of a lipid,
            return its histogram with edges."""
            x = array[:, 0]
            y = array[:, 1]

            lipid, e1, e2 = np.histogram2d(x, y, density=True, bins=bins)

            return lipid, e1, e2

        # Plot histogram image
        H, xe, ye = lipid_density(l_xyz, number.value)
        minx = np.abs(xe.min())
        miny = np.abs(ye.min())

        p1 = figure(plot_height=640, plot_width=640,  tools='')
        image_source = ColumnDataSource(data=dict(image=[]))
        img = p1.image(image="image", x=xe[0], y=ye[0], dw=xe[-1] + minx, dh=ye[-1] + miny, palette=bokehpalette, source=image_source)

        circle_prot_source = ColumnDataSource(data=dict(x1=prot_xyz[:, 0], y1=prot_xyz[:, 1]))
        p1.circle(x="x1", y="y1", source=circle_prot_source, size=2, fill_alpha=0.2)

        cb_palette = LinearColorMapper(palette=bokehpalette, low=H.min(), high=H.max())
        color_bar = ColorBar(color_mapper=cb_palette, width=8,  location=(0,0), label_standoff=10)
        color_bar.formatter = BasicTickFormatter(use_scientific=False)
        p1.add_layout(color_bar, 'right')

        # Make graph pretty
        p1.xgrid.grid_line_color = None
        p1.ygrid.grid_line_color = None
        p1.xaxis.major_tick_line_color = None
        p1.xaxis.minor_tick_line_color = None
        p1.yaxis.major_tick_line_color = None
        p1.yaxis.minor_tick_line_color = None
        p1.xaxis.major_label_text_font_size = '0pt'
        p1.yaxis.major_label_text_font_size = '0pt'
        p1.grid.visible = False
        p1.toolbar.logo = None
        p1.toolbar_location = None


        def update_all(cond=False, cmap=False):
            """
            Update the image showing all proteins.
            """
            if cond:
                # For efficiency execute only if GPCR structure changes
                if filenames is not None:
                    update_l_xyz = np.load(lipid.value)
                else:
                    update_l_xyz = flat_density(t, [lipid.value], t.n_frames)

                x_min, x_max = update_l_xyz[:, 0].min(), update_l_xyz[:, 0].max()
                y_min, y_max = update_l_xyz[:, 1].min(), update_l_xyz[:, 1].max()
                z_min, z_max = update_l_xyz[:, 2].min(), update_l_xyz[:, 2].max()
                zrange.start = z_min
                zrange.end   = z_max

                index = np.where((update_l_xyz[:, 2] > zrange.value[0]) &
                                (update_l_xyz[:, 2] < zrange.value[1]))
                l_xyz_new = update_l_xyz[index]
            else:
                l_xyz_new = l_xyz

            if cmap:
                # For efficiency execute only if image colormap changes
                cb_cut_value = 256 - cbar_range.value

                cmap = color_map.value
                colormap = cm.get_cmap(cmap)
                bokehpalette = [mpl.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
                bp_i = 0
                while bp_i < len(bokehpalette[:cb_cut_value]):
                    bokehpalette[bp_i] = '#ffffff'
                    bp_i += 1

                img.glyph.color_mapper.palette = bokehpalette
                color_bar.color_mapper.palette = bokehpalette

            # Update histogram image
            H, xe, ye = lipid_density(l_xyz_new, number.value)
            minx = np.abs(xe.min())
            miny = np.abs(ye.min())

            img.glyph.dw = xe[-1] + minx
            img.glyph.dh = ye[-1] + miny

            # update image source
            image_source.data = dict(image=[H])

        def update_protein():
            """
            Update the protein representation.
            """

            if protein.value == "Yes":
                circle_prot_source.data = dict(
                    x1=prot_xyz[:, 0],
                    y1=prot_xyz[:, 1]
                )

            elif protein.value == "No":
                circle_prot_source.data = dict(
                    x1=[],
                    y1=[]
                )

        def update_cbar():

            cb_cut_value = 256-cbar_range.value

            cmap = color_map.value
            colormap = cm.get_cmap(cmap)
            bokehpalette = [mpl.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]

            bp_i = 0
            while bp_i < len(bokehpalette[:cb_cut_value]):
                bokehpalette[bp_i] = '#ffffff'
                bp_i += 1

            img.glyph.color_mapper.palette = bokehpalette
            color_bar.color_mapper.palette = bokehpalette


        # event listeners
        controls = [lipid, zrange]
        for control in controls:
            control.on_change('value', lambda attr, old, new: update_denstype(cond=True))

        number.on_change('value', lambda attr, old, new: update_denstype())
        color_map.on_change('value', lambda attr, old, new: update_denstype(cmap=True))
        protein.on_change('value', lambda attr, old, new: update_protein())
        cbar_range.on_change('value', lambda attr, old, new: update_cbar())

        # deal with what gets updated and what not.
        def update_denstype(cond=False, cmap=False):
            update_all(cond, cmap)
            update_protein()

        update_denstype()
        input1 = row([lipid, color_map, protein])
        input2 = row([p1, cbar_range])
        input3 = row([number, zrange])
        input3 = column([input1, input2, input3])

        l = layout([input3])

        doc.add_root(l)
        doc.title = "Density App"

    return show(densityApp)
