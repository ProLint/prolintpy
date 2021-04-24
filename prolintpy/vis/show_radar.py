import yaml
import numpy as np
import pandas as pd

# bokeh aplication
from bokeh.document import Document

from bokeh.plotting import figure
from bokeh.layouts import layout, column
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.models.widgets import Select
from bokeh.palettes import RdBu

from bokeh.io import show, output_notebook
from bokeh.themes import Theme

from prolintpy.utils.misc import unit_poly_verts, radar_patch, star_curv

def show_radar(df=None, filename=None):
    """Visualize lipid-protein contacts using a Bokeh radar map.

    This is useful if you want to (i) compare the binding profile of two residues and
    (ii) when you want to see how different contact metrics compare against each other.

    Parameters
    ----------

    df : pd.DataFrame
         A contacts df. This will be the output of pl.contacts_dataframe() command.

    filename: string
        A string pointing to the location of a csv file containing data that can be loaded into a pandas DataFrame.

    """
    if df is None and not filename:
        raise Exception("show_radar can take either a pandas DataFrame or a csv file that can be loaded into a pandas DataFrame.")
    elif filename:
        df = pd.read_csv(filename)

    output_notebook()

    def radarApp(doc):

        nonlocal df
        df = df.copy()
        for col in df.columns:
            if col.endswith('Error'):
                del df[col]

        gpcrs = list(df.Protein.unique())
        lipids = list(df.Lipids.unique())
        radius = list(df.Radius.unique())
        radius = [str(r) for r in radius]

        # Normalize [0-1] all contact parameters with respect to each other.
        for cmn in list(df.columns)[:-5]:
            for lipid in lipids:
                for radius_value in radius:
                    radius_value = float(radius_value)

                    dfmin = df[
                        (df.Lipids == lipid) &
                        (df.Radius == radius_value)
                        ][cmn].min()

                    dfmax = df[
                        (df.Lipids == lipid) &
                        (df.Radius == radius_value)
                        ][cmn].max()

                    lip_index = list(df[
                        (df.Lipids == lipid) &
                        (df.Radius == radius_value)
                        ].index)

                    df.loc[lip_index, cmn] = df[
                        (df.Lipids == lipid) &
                        (df.Radius == radius_value)
                        ][cmn].apply(lambda x: (x-dfmin) / (dfmax - dfmin))

        df = df.fillna(0)


        # 1 unique set of residues - works also if there's only 1 set available.
        resname = df[(df.Protein == gpcrs[0]) & (df.Lipids == lipids[0]) & (df.Radius == float(radius[0]))].ResName.to_list()
        resid   = df[(df.Protein == gpcrs[0]) & (df.Lipids == lipids[0]) & (df.Radius == float(radius[0]))].ResID.to_list()

        residues = []
        for rn, ri in zip(resname, resid):
            residues.append("{}-{}".format(rn, ri))


        import random
        # Create Input controls
        gpcr1 = Select(title="Protein", value=gpcrs[0], options=gpcrs, width=100)
        gpcr2 = Select(title="Protein", value=gpcrs[0], options=gpcrs, width=100)
        lipid1 = Select(title="Lipid", value=lipids[0], options=lipids, width=100)
        lipid2 = Select(title="Lipid", value=lipids[0], options=lipids, width=100)
        residue1 = Select(title="Residue", value=random.choice(residues), options=residues, width=120)
        residue2 = Select(title="Residue", value=random.choice(residues), options=residues, width=120)
        radius1 = Select(title="Radius", value=radius[0], options=radius, width=100)
        radius2 = Select(title="Radius", value=radius[0], options=radius, width=100)


        # num_vars = 6
        # cols_to_display = [x for x in list(df.columns) if not x.endswith('Error')]
        # num_vars = len(cols_to_display)

        num_vars = len(list(df.columns)[:-5])
        centre = 0.5

        theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)
        theta += np.pi/2

        verts = unit_poly_verts(theta, centre)
        xv = [v[0] for v in verts]
        yv = [v[1] for v in verts]

        f1 = df[
                (df.Protein == str(gpcr1.value)) &
                (df.Lipids == str(lipid1.value)) &
                (df.ResID == int(str(residue1.value).split("-")[1])) &
                (df.Radius == float(radius1.value))
                ].loc[:, list(df.columns)[:num_vars]].to_numpy().reshape(num_vars,)

        f2 = df[
                (df.Protein == str(gpcr2.value)) &
                (df.Lipids == str(lipid2.value)) &
                (df.ResID == int(str(residue2.value).split("-")[1])) &
                (df.Radius == float(radius2.value))
                ].loc[:, list(df.columns)[:num_vars]].to_numpy().reshape(num_vars,)

        flist = [f1, f2]

        p = figure(plot_height=400, plot_width=800, x_range=(-0.3, 1.3), y_range=(-0.3, 1.3), tools="")
        cmap = RdBu[num_vars]
        colors = cmap[:1] + cmap[-1:]

        patch1 = ColumnDataSource(data=dict(xi=[], yi=[]))
        circle1 = ColumnDataSource(data=dict(xt=[], yt=[]))
        xt, yt = radar_patch(f1, theta, centre)
        xi, yi = star_curv(xt + xt[:1], yt + yt[:1])
        patch = p.patch(x="xi", y="yi", fill_alpha=0.15, fill_color=colors[0], line_color=colors[0], source=patch1)
        circle = p.circle(x='xt', y='yt', size=10, alpha=0.5, line_color=colors[0], fill_color=colors[0], source=circle1)

        patch2 = ColumnDataSource(data=dict(xi=[], yi=[]))
        circle2 = ColumnDataSource(data=dict(xt=[], yt=[]))
        xt, yt = radar_patch(f2, theta, centre)
        xi, yi = star_curv(xt + xt[:1], yt + yt[:1])
        patch = p.patch(x="xi", y="yi", fill_alpha=0.15, fill_color=colors[1], line_color=colors[1],  source=patch2)
        circle = p.circle(x='xt', y='yt', size=10, alpha=0.5, line_color=colors[1], fill_color=colors[1], source=circle2)

        # Draw grid so the plot looks like a polar graph.
        text_label = list(df.columns)[:num_vars] + ['']
        p.circle(x=0.5, y=0.5, size=20, fill_alpha=0, line_color='black', line_dash='dashed', line_alpha=0.2)
        for size in [x*100 for x in range(1, 10)]:
            p.circle(x=0.5, y=0.5, size=size, fill_alpha=0, line_color='black', line_dash='dashed', line_alpha=0.2)

        line = np.linspace(-2, 2, 100)
        p.line(line, (line * 0.58)+0.21, line_color='black', line_dash='dashed', line_alpha=0.2)
        p.line(line, (line * (-0.58))+0.79, line_color='black', line_dash='dashed', line_alpha=0.2)
        p.line([0.5]*100, line, line_color='black', line_dash='dashed', line_alpha=0.2)

        # Hide axes.
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        p.xaxis.major_tick_line_color = None
        p.xaxis.minor_tick_line_color = None
        p.yaxis.major_tick_line_color = None
        p.yaxis.minor_tick_line_color = None
        p.xaxis.major_label_text_font_size = '0pt'
        p.yaxis.major_label_text_font_size = '0pt'

        p.axis.visible = None
        p.toolbar.logo = None
        p.toolbar_location = None

        # Draw the hexagon with labels attached.
        p.line(x=xv+[centre], y=yv+[1], color='black')

        xv[1] = xv[1] - 0.13
        xv[2] = xv[2] - 0.1
        xv[-2] = xv[-2] + 0.16
        xv[-1] = xv[-1] + 0.1
        # yv[3] = yv[3] - 0.06

        label_source = ColumnDataSource({'x':xv + [centre ],'y':yv + [1],'text':text_label})
        labels = LabelSet(x="x",y="y",text="text",source=label_source, level='glyph', text_align="center")
        p.add_layout(labels)

        def select_dfGPCR():
            gpcr1_val, gpcr2_val  = gpcr1.value, gpcr2.value
            lipid1_val, lipid2_val = lipid1.value, lipid2.value
            res1_val, res2_val = str(residue1.value).split("-")[1], str(residue2.value).split("-")[1]

            f1 = df[
                (df.Protein == str(gpcr1_val)) &
                (df.Lipids == str(lipid1_val)) &
                (df.ResID == int(res1_val)) &
                (df.Radius == float(radius1.value))
                ].loc[:, list(df.columns)[:num_vars]].to_numpy().reshape(num_vars,)
            f2 = df[
                (df.Protein == str(gpcr2_val)) &
                (df.Lipids == str(lipid2_val)) &
                (df.ResID == int(res2_val)) &
                (df.Radius == float(radius2.value))
                ].loc[:, list(df.columns)[:num_vars]].to_numpy().reshape(num_vars,)

            return (f1, f2)

        def update1():
            gpcr_val = gpcr1.value
            lipid1_val = lipid1.value
            radius1_val = radius1.value

            resname = df[(df.Protein == str(gpcr_val)) & (df.Lipids == str(lipid1_val)) & (df.Radius == float(radius1_val))].ResName.to_list()
            resid   = df[(df.Protein == str(gpcr_val)) & (df.Lipids == str(lipid1_val)) & (df.Radius == float(radius1_val))].ResID.to_list()

            residues = []
            for rn, ri in zip(resname, resid):
                residues.append("{}-{}".format(rn, ri))

            residue1.options = residues
            if residue1.value not in residues:
                residue1.value = residues[0]

            f1 = select_dfGPCR()[0]
            xt, yt = radar_patch(f1, theta, centre)
            xi, yi = star_curv(xt + xt[:1], yt + yt[:1])

            circle1.data=dict(xt=xt, yt=yt)
            patch1.data=dict(xi=xi, yi=yi)

        def update2():
            gpcr_val = gpcr1.value
            lipid2_val = lipid2.value
            radius2_val = radius2.value

            resname = df[(df.Protein == str(gpcr_val)) & (df.Lipids == str(lipid2_val)) & (df.Radius == float(radius2_val))].ResName.to_list()
            resid   = df[(df.Protein == str(gpcr_val)) & (df.Lipids == str(lipid2_val)) & (df.Radius == float(radius2_val))].ResID.to_list()

            residues = []
            for rn, ri in zip(resname, resid):
                residues.append("{}-{}".format(rn, ri))

            # options = [str(x) for x in options]
            residue2.options = residues
            if residue2.value not in residues:
                residue2.value = residues[0]

            f2 = select_dfGPCR()[1]
            xt, yt = radar_patch(f2, theta, centre)
            xi, yi = star_curv(xt + xt[:1], yt + yt[:1])

            circle2.data=dict(xt=xt, yt=yt)
            patch2.data=dict(xi=xi, yi=yi)


        controls1 = [gpcr1, lipid1, residue1, radius1]
        controls2 = [gpcr2, lipid2, residue2, radius2]

        for control1 in controls1:
            control1.on_change('value', lambda attr, old, new: update1())
        for control2 in controls2:
            control2.on_change('value', lambda attr, old, new: update2())

        inputs1 = column([gpcr1, lipid1, residue1, radius1])
        inputs2 = column([gpcr2, lipid2, residue2, radius2])

        l = layout([[inputs1, p, inputs2]])

        update1()
        update2()

        doc.add_root(l)

        doc.title = "Radar App"
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

    return show(radarApp)
