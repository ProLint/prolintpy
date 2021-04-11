import yaml
import numpy as np

from bokeh.document import Document

from bokeh.plotting import figure
from bokeh.layouts import layout
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Select
from bokeh.layouts import row
from bokeh.io import show, output_notebook

from bokeh.themes import Theme

from bokeh.models import Range1d


def show_distances(distances, norm_factor=1000000):
    """Visualize the distance between protein residues and lipids.

    Parameters
    ----------
    distances : tuple
        This is the output of a pl.ComputeContacts.compute_distances() call.
        See docs for more help.

    norm_factor : int
        Normalize time dimension (i.e. convert time to nano or microsecond).
        Default is 1000000 and will convert picoseconds to microseconds.

    """

    output_notebook()

    def distanceApp(doc):

        global df

        residues = list(distances.keys())
        df, meta = distances[residues[0]]
        residues = [str(x) for x in residues]

        time = df['Time'].to_list()

        proteins = [meta['protein']]
        protein_id = meta['protein_id']
        protein_id = [str(x) for x in protein_id]
        min_val = meta['min_val']

        columns = []
        for resid in residues:
            for col in distances[int(resid)][0].columns[1:]:
                columns.append(col)

        source_dict = {}
        for i in np.unique(columns):
            source_dict[i] = []
        source_dict['x'] = []

        source = ColumnDataSource(data=source_dict)

        protein_name = Select(title="Protein", value=proteins[0], options=proteins, width=100)
        pc = Select(title="protein_id", value=protein_id[0], options=protein_id, width=100)
        res = Select(title="Residue Selection", value=residues[0], options=residues, width=100)

        p = figure(plot_width=1500, plot_height=400, tools='pan, box_zoom, ywheel_zoom, save, reset, help')

        for i in np.unique(columns):
            color = '#%02x%02x%02x' % tuple(np.random.choice(range(256), size=3))
            p.line(x='x', y=i, line_color=color, source=source, line_width=2)

        p.y_range = Range1d(min_val-(min_val*0.15), 2.5)
        p.x_range = Range1d(time[0], time[-1])

        p.toolbar.autohide = True
        p.axis.axis_label_text_font_size = "12pt"
        p.axis.axis_label_text_font_style = "bold"
        p.title.align = 'center'

        p.xaxis.axis_label = "Trajectory Time"
        p.yaxis.axis_label = "Distance (nm)"


        def update():

            df = distances[int(res.value)][0]
            time = df['Time'].to_list()
            df['NaNs'] = ['NaN'] * len(time)

            show_columns = [x for x in df.columns if x.startswith(pc.value + '_')]
            not_columns = [x for x in np.unique(columns) if x not in show_columns]

            new_source_dict = {}
            new_source_dict['x'] = df['Time']
            for col in show_columns:
                new_source_dict[col] = df[col]

            for col in not_columns:
                new_source_dict[col] = df['NaNs']

            source.data = new_source_dict


        controls = [protein_name, pc, res]
        for control in controls:
            control.on_change('value', lambda attr, old, new: update())

        inputs2 = row([*controls])
        layout1 = layout([[inputs2]])
        layout2 = layout([p])

        update()

        doc.add_root(layout1)
        doc.add_root(layout2)

        doc.title = "Distance Calculations"
        doc.theme = Theme(json=yaml.load("""
            attrs:
                Figure:
                    toolbar_location: above
                    height: 500
                    width: 800
                Grid:
                    grid_line_dash: [6, 4]
                    grid_line_color: black
        """, Loader=yaml.FullLoader))

    return show(distanceApp)
