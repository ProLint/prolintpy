# prolintpy
prolintpy is a lightweight python library that aimes to automate the analysis and visualization of **Pro**tein-**L**ipid **int**eractions.

The main goal and purpose of all the tools that are distributed as part of the <a href="https://prolint.ca" target="_blank">ProLint</a> framework is
to bridge this widening gap between **data generation** and **gaining insight** on biologically-relevant interactions between lipids and proteins.

prolintpy is the library that the ProLint webserver uses on the backend to automate topology generation and analysis of lipid-protein interactions. Please note, however, that
prolintpy includes a dedicated interface for the visualization of lipid-protein interactions similar to the webserver.

## What does this tool do?
TLDR: you can use `prolintpy` for the following:
<ol>
<li>Automatically generate a topology description of your system (no <span style="font-style: oblique;">tpr</span> file needed)</li>
<li>Calculate contact-based metrics for lipid-protein interactions</li>
<li>Calculate 2D and 3D densities (in progress)</li>
<li>Calculate physics-based properties (in progress)</li>
<li>Interactively visualize lipid-protein interactions</li>
</ol>


## Installation

To install `prolintpy` simply execute: 

```sh
python -m pip install prolintpy
```

This should work on most systems.
On Windows and even WSL 1, MDTraj may present a problem to install. In that case, you may want to use `conda` to 
install MDTraj first: 

```sh
# create new environment
conda create --name prolint python=3.7
conda install -c conda-forge mdtraj
python -m pip install prolintpy
```

### Installing from source

If you want to install directly from the github repository then you can do that by typing: 

```sh
git clone https://github.com/ProLint/prolintpy.git
cd prolintpy
python setup.py install
```

If you are using Windows, the same thing mentioned above applies. 

## Getting Started

Please follow the instructions provided in the <a href="https://prolint.github.io/prolintpy" target="_blank">documentation</a> to get started. Note that, to use the visualization interface of prolintpy,
you should use JupyterLab. At the top of your notebook file, make sure to call the `output_notebook` function:

```python
from bokeh.io import output_notebook
output_notebook()
```

Additionally, if you want to use the `show_contact_projection` function, make sure that your installation of `nglview` is working properly.
Follow the instruction provided <a href="https://github.com/nglviewer/nglview" target="_blank">there</a> to ensure your installation is running correclty.

## Data Files
Before you load the data to prolintpy make sure to first remove water & ions beads. Leave only membrane and protein beads in the system.

