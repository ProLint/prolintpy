import yaml
import numpy as np
import pandas as pd
import colorcet as cc

from bokeh.document import Document

from bokeh.plotting import figure
from bokeh.layouts import layout
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, Select, TextInput
from bokeh.layouts import row
from bokeh.palettes import Reds
from bokeh.transform import linear_cmap
from bokeh.io import show, output_notebook

from bokeh.themes import Theme

def show_points(df=None, filename=None, **kwargs):
    """Visualize lipid-protein contacts using a Bokeh scatter plot.

    Parameters
    ----------

    df : pandas.DataFrame
        A contacts dataframe. This will be the output of pl.contacts_dataframe() command.

    filename: string
        A string pointing to the location of a csv file containing data that can be loaded into a pandas DataFrame.

    size : int
        The size of points to dispaly. Default is 10.

    """
    if df is None and not filename:
        raise Exception("show_points can take either a pandas DataFrame or a csv file that can be loaded into a pandas DataFrame.")
    elif filename:
        df = pd.read_csv(filename)

    output_notebook()

    def pointsApp(doc):

        radii = [str(x) for x in df.Radius.unique()]

        # Create Input controls
        number = Slider(title="Value Cutoff", value=0, start=0, end=4, step=0.1, width=150)

        residue = TextInput(title="Residue name (3 letter code):", width=200)

        gpcr = Select(title="Proteins", value=list(df.Protein.unique())[0],
                    options=list(df.Protein.unique()), width=100)

        lipid = Select(title="Lipids", value=list(df.Lipids.unique())[0],
                    options=list(df.Lipids.unique()), width=100)

        radius = Select(title="Radius", value=radii[-1], options=radii, width=100)

        options = list(df.columns)[:-5] + ['ResID']
        x_axis = Select(title="X Axis", options=options, value="ResID", width=150)
        y_axis = Select(title="Y Axis", options=options, value=options[0], width=150)


        cc_colors = [x for x in cc.all_original_names() if x.startswith('linear') or x.startswith('rainbow')]
        cmap = Select(title="Colormap", options=cc_colors, value='linear_kryw_0_100_c71', width=150)

        # Create Column Data Source that will be used by the plot
        source = ColumnDataSource(data=dict(x=[], y=[], ResName=[], ResID=[], Protein=[]))

        TOOLTIPS=[
            ("ResName", "@ResName"),
            ("ResID", "@ResID"),
            ("Value", "@y")
        ]

        mapper = linear_cmap(field_name='y', palette=cc.CET_L19,
                            low=df[df.Protein == gpcr.value][y_axis.value].min(),
                            high=df[df.Protein == gpcr.value][y_axis.value].max())

        p = figure(tooltips=TOOLTIPS,)

        global c
        c = p.circle(x="x", y="y", source=source, line_color='black', fill_color=mapper, **kwargs)

        p.toolbar.autohide = True
        p.axis.axis_label_text_font_size = "12pt"
        p.axis.axis_label_text_font_style = "bold"
        p.title.align = 'center'

        def update(df):
            y_value = y_axis.value
            x_value = x_axis.value

            df = df[
                (df[y_value] >= number.value) &
                (df['Protein'] == gpcr.value) &
                (df['Lipids'] == lipid.value) &
                (df['Radius'] == float(radius.value))
            ]
            if (residue.value != ""):
                df = df[df.ResName.str.contains(residue.value.upper())==True]

            mapper = linear_cmap(field_name='y', palette=cc.palette[cmap.value],
                                low=df[df.Protein == gpcr.value][y_value].min(),
                                high=df[df.Protein == gpcr.value][y_value].max())

            c.glyph.fill_color = mapper

            p.xaxis.axis_label = x_value
            p.yaxis.axis_label = y_value
            p.title.text = "Showing %d Data Points  " % len(df)

            source.data = dict(
                x=df[x_value],
                y=df[y_value],
                ResName=df["ResName"],
                ResID=df["ResID"],
                Protein=df["Protein"],
            )

        controls = [number, gpcr, lipid, radius, y_axis, x_axis, residue, cmap]
        for control in controls:
            control.on_change('value', lambda attr, old, new: update(df))

        sizing_mode = 'scale_width'

        inputs = row(*controls, sizing_mode=sizing_mode)
        inputs2 = row([gpcr, lipid, radius, residue], sizing_mode=sizing_mode)
        inputs3 = row([number, x_axis, y_axis, cmap], sizing_mode=sizing_mode)

        layout1 = layout([[inputs2]], sizing_mode=sizing_mode)
        layout2 = layout([p], sizing_mode=sizing_mode)
        layout3 = layout([inputs3], sizing_mode="scale_width")

        update(df)

        doc.add_root(layout1)
        doc.add_root(layout2)
        doc.add_root(layout3)
        doc.title = "Scatter Application"
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

    return show(pointsApp)