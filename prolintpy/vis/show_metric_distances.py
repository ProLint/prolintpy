import yaml
import numpy as np
import pickle

from bokeh.document import Document

from bokeh.plotting import figure
from bokeh.layouts import layout
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Select
from bokeh.layouts import row
from bokeh.io import show, output_notebook

from bokeh.models import Range1d
from bokeh.themes import Theme


def show_metric_distances(lipid_residue_distances):
    """Visualize the distance between protein residues and lipids.

    """

    output_notebook()

    def distanceApp(doc: Document) -> None:

        distances = lipid_residue_distances[0]
        time = lipid_residue_distances[1]['time']
        proteins = lipid_residue_distances[1]['protein']
        lipids = list(distances.keys())
        residues = list(distances[lipids[0]].keys())
        protein_id = list(distances[lipids[0]][residues[0]].keys())

        residues = [str(x) for x in residues]
        protein_id = [str(x) for x in protein_id]

        source = ColumnDataSource(data=dict(x=[], y=[]))

        protein_name = Select(title="Protein", value=proteins[0], options=proteins, width=100)
        pc = Select(title="Protein Copy", value=protein_id[0], options=protein_id, width=100)
        res = Select(title="Residue Selection", value=residues[0], options=residues, width=100)
        lip = Select(title="Lipid Selection", value=lipids[0], options=lipids, width=100)
        p = figure(plot_width=1500, plot_height=400, tools='pan, box_zoom, ywheel_zoom, save, reset, help')

        color = '#%02x%02x%02x' % tuple(np.random.choice(range(256), size=3))
        p.line(x='x', y='y', line_color=color, source=source, line_width=4)

        p.y_range = Range1d(0, 4)
        p.x_range = Range1d(time[0], time[-1])

        p.toolbar.autohide = True
        p.axis.axis_label_text_font_size = "12pt"
        p.axis.axis_label_text_font_style = "bold"
        p.title.align = 'center'

        p.xaxis.axis_label = "Trajectory Time"
        p.yaxis.axis_label = "Distance (nm)"

        def update():

            protein_selected = protein_name.value
            copy_selected = pc.value
            residue_selected = res.value
            lipid_selected = lip.value

            residues = list(distances[lipid_selected].keys())
            if int(residue_selected) not in residues:
                residue_selected = residues[0]

            res.options = [str(x) for x in residues]

            y = distances[lipid_selected][int(residue_selected)][int(copy_selected)]
            source.data = dict(
                x=time,
                y=y
            )

        controls = [protein_name, pc, res, lip]
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
                    width: 1000
                Grid:
                    grid_line_dash: [6, 4]
                    grid_line_color: black
        """, Loader=yaml.FullLoader))

    return show(distanceApp)

