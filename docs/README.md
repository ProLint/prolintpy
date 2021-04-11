# prolintpy
prolintpy is a lightweight python library that aimes to automate the analysis and visualization of **Pro**tein-**L**ipid **int**eractions.

The main goal and purpose of all the tools that are distributed as part of the <a href="https://prolint.readthedocs.io" target="_blank">ProLint</a> framework is
to bridge this widening gap between **data generation** and **gaining insight** on biologically-relevant interactions between lipids and proteins.

prolintpy is the library that the ProLint webserver uses on the backend to automate topology generation and analysis of lipid-protein interactions. Please note, however, that
prolintpy includes a dedicated interface for the visualization of lipid-protein interactions similar to the webserver.

## What does this tool do?
TLDR: you can use `prolintpy` for the following:
<ol>
<li>Automatically generate a topology description of your system (no <span style="font-style: oblique;">tpr</span> file needed)</li>
<li>Calculate contact-based metrics for lipid-protein interactions</li>
<li>Calculate 2D and 3D densities</li>
<li>Calculate physics-based properties</li>
<li>Interactively visualize lipid-protein interactions</li>
</ol>


## Installation

Getting `prolintpy` is simple, especially if you have `conda` installed.

```sh
# create new environment
conda create --name prolint python=3.7
conda install -c conda-forge mdtraj
```
and then install `prolintpy` using `pip`:
```sh
python -m pip install prolint
```

### Installing from source
If you want to install directly from the github repository then you can do that by typing:

```sh
# create new environment
conda create --name prolint python=3.7
# install dependencies
conda install -c conda-forge numpy pandas==0.24.0 mdtraj scipy pyyaml colorcet bokeh==1.4.0 networkx nglview==2.7.7 matplotlib jupyterlab
```

After that, you clone this directory and install it, using:

```sh
git clone https://github.com/bisejdiu/prolint.git
cd prolint
python setup.py install
```

## Getting Started

Please follow the instructions provided in the <a href="https://prolint.readthedocs.io" target="_blank">documentation</a> to get started. Note that, to use the visualization interface of prolintpy,
you should use JupyterLab. At the top of your notebook file, make sure to call the `output_notebook` function:

```python
from bokeh.io import output_notebook
output_notebook()
```

Additionally, if you want to use the `show_contact_projection` function, make sure that your installation of `nglview` is working properly.
Follow the instruction provided <a href="https://github.com/nglviewer/nglview" target="_blank">there</a> to ensure your installation is running correclty.

## Data Files
Before you load the data to prolintpy make sure to first remove water & ions beads. Leave only membrane and protein beads in the system.

