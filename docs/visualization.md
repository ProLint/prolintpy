# Visualization Reference

`prolintpy` automates analysis and visualization of lipid-protein interactions. <br>
This section outlines how to use ProLint for for visualizing contact information using interactive applications.

In contrast to analysis, data visualization is more variable and different users will
have different preferences for how to best display results. `prolintpy` does not make any assumptions about the data or the preferences of the user.

Instead, the main objective is to allow users to quickly and without hassle obtain valuable insight on
protein-lipid interaction details in their simulated systems. To do this, it uses dedicated, modern,
popular and extremely powerful libraries that do that heavy
lifting in the background. All visualization tools supported by `prolintpy` are interactive and
are intended to give users complete access to their data.

Make sure you are using the JupyterLab environmet for this, and at the beginning execute:

```python
from bokeh.io import output_notebook
output_notebook()
```

### Point distribution

Once you have build a pandas DataFrame from the calculated contacts (by, for instance, running the `contacts_dataframe()` function),
you can provide it as input to the appropriate `prolintpy` visualization apps:

```python
pl.show_points(df, size=15)
```

When runing in a JupyterLab environment, will give the following output application:

<img src="assets/images/points.png"></img>

All of the provided tools (sliders, text boxes and select boxes) are intuitive and easy to understand. Importantly, the data that you
visualize does not strictly have to have been calculated using prolintpy and you can use any other contacts you like. The app requires the presence
of these columns: ResID, ResName, Protein, Lipids, Radius put at the end of the DataFrame.
All are strings, except Radius which must contain be float values. If you have a different dataframe that you wish to visualize,
either rename your columns to match the required names, or simply add placeholder value if you do not require them.

### Metric comparison

Because by default ProLint calculates several contact metrics and users may add custom definitions,
you may want to see how these parameters compare against each other. In particular, this becomes
important when you have two or more residues that interact strongly with a lipid, and you want to
compare their interactions against each other. You can do this using ProLint by typing::

```python
pl.show_radar(df)
```

This will normalize all contact metrics (there should be at least three present) and use a radar plot
for the visualization:

<img src="assets/images/radar.png"></img>

### Distance calculations 1

A very important (perhaps the most important?) calculation that is commonly done in lipid-protein interaction studies is measuring the
distance between a residue and a lipid as a function of simulation time. This gives you a clear idea if the
lipid is interacting preferentially with a residue or not. `prolintpy` provides two different ways to get distance information on
lipid-protein interactions. The first method, presented in this section, is automated and relies on the prior calculation of contact-based metrics.

The way it works is that it goes through the calculated metrics, sorts them, and gets the top-ranking residues and lipids. It then goes over each
residue and lipid combination and gets the best/strongest contact (that is, the contact that is maintained most strongly between the specific residue and lipid).

Here is an example application:

```python
from prolintpy.core.computecontacts import retrieve_distances

# Note that the only calculated value needed is the contacts dataframe.
# The other arguments are simply definitions about your system.
distances_dict, SYSTEM_LIPIDS, lipids_found = retrieve_distances(protein_dataframe, group_lipids=False, resolution="martini", lipids=lipids, top_nr=30)
distances = contacts.compute_lipid_distances(t, proteins[0], distances_dict, SYSTEM_LIPIDS, lipids_found)

# Visualize results
pl.show_metric_distances(distances)
```

Note that the `top_nr` argument specifies the number of top ranking results based on contact-metrics to consider. You can increase it to get a larger pool of interactions.<br>
This will visualize contacts found using the following application:

<img src="assets/images/mdistances.png"></img>

Note how you can change the Lipid Selected and this will update the list in the Residue Selection dropdown. Similarly, if you have multiple proteins in your
system, the values displayed will be grouped acording to the protein copy. Note that if you did not merge proteins, you need to undo that. See here:
https://github.com/ProLint/ProLint/blob/main/prolint/calcul/tasks.py#L141 for more information on how this is done.

### Distance calculations 2

The second way `prolintpy` calculates and visualizes distances is by not relying on any prior calcualted metrics. Instead, you simply supply the protein
and list of residues along with threshold arguments, and `prolintpy` will then calculate distance measurements.

Given a list of input residues,
this function will loop through all the lipids in the system and display distances with best ranking lipids. Ranking is decided
based on the following parameters:

| Argument      | Default | Description                                                                                                              |
| ------------- | ------- | ------------------------------------------------------------------------------------------------------------------------ |
| distance_co   | 0.7     | A cutoff distance (nm) that a lipid must satisfy for `percentile_co` frames of the trajectory.                           |
| percentile_co | 0.05    | The percentage of the trajectory (measured in frames) that a lipid must be within the `distance_co` for it to be stored. |

Here is an example application:

```python
# Calculate contacts for the first protein in the system and residues in the range 80-119 with cholesterol
dist = contacts.compute_distances(t, proteins[0], [*range(1, 10)], 'CHOL', 'ROH', percentile_co=0.05, distance_co=0.7)
pl.show_distances(dist)
```

You'll get an output like the following (note that the actual application also has a dropdown menu to select the residue which is not shown in this image):

<img src="assets/images/distance.png"></img>

Note how you may and usually will, get multiple lipids satisfying the threshold parameters. This is a really cool way to visualize interactions along
long duration of the trajectory and observe many lipid binding and unbinding events.

### Heatmap & Density Viewer

NGL Viewer based application to visualize contact metrics as contact heatmaps that are projected on the
surface of the protein. You visualize it in `prolintpy` using:

```python
# depending on the size of the protein this may take a couple of seconds to render.
# if you have multiple cutoffs then you also need to filter the dataframe using one value
contact_values = df[df.Lipids == "CHOL"].Longest_Duration.to_list()
pl.show_contact_projection(t, bf=contact_values, protein=proteins[0], cmap="Reds")
```

In the first argument we give information about the topology of the system, in the second we provide
the contacts and in the third argument we tell `prolintpy` which protein to use. The key here is that the
size of `bf` has to match the number of residues in `protein`.
We can also use the `ngl_repr` argument to specify the structure representation. For example, if we use options like:
`surface`, `cartoon`, or `ball+stick` we get the following outputs:

<img src="assets/images/projection.png"></img>

Please note that `prolintpy` uses the nglview library. So this application is intendent to cover a majority
of use cases when it comes to visualizing contacts with lipids, but obviously it is not intendent to be
comprehensive or provide many options.

Also note that currently only `martini` simulations are supported for this visualization and support for `atomistic` input files
will be added very soon.

### Graph Network

Lipids and proteins are represented as nodes and their interactions with each other as edges.
The width of the latter corresponds to the degree of the underlying interaction. To visualize
lipid-protein interactions as a graph network you execute the following commands:

```python
# df.columns[2] indicates which values we want to visualize
pl.network(lipids, df, df.columns[2], grouped=True)
pl.show_network('Protein0_Longest_Duration_0.5_network.json')  # may have to execute twice.
```

In the first command, we store the interactions in a json file and then visualize them in the second
command. This will give an output like this:

<img src="assets/images/network.png"></img>

The size of lipid nodes corresponds to the ratio of that particular lipid group in simulated systems.
The size of residue nodes corresponds to the relative content of each residue in the protein being shown.

### Density Maps

The positions visited during the simulation by lipids are bined using a 2D histogram and colored using a color map.
This shows the preferential localization of lipids. This type of analysis depends on the protein being centered, which
may have been fixed as part of the system setup or can be centered using other trajectory processing tools.
ProLint allows you to calculate 2D densities using a very simple syntax:

```python
# Show the density of cholesterol and POPS lipids for a martini simulation:
pl.show_density(t, ['CHOL', 'POPS'], "martini")
```

The command above will yieald a result like the following:

<img src="assets/images/density_map.png"></img>

There is a lipid selectbox to change the lipid being displayed and a colormap as well as Show Protein select box to change the
type of color gradient used to visualize densities and show/hide the protein respectively. There are three sliders: to change the
number of bins used when binning the data (the more bins the finer the resolution), a z-axis range that allows you to specify which
range of the bilayer to use for the density calculation (e.g. only the upper leaflet), and a colorbar slider to hide values under
a threshold.
