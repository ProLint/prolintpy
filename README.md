# prolintpy
A python package for the automated analysis and visualization of **Pro**tein-**L**ipid **int**eractions.


prolintpy is a lightweight python library that is used by the ProLint webserver on the backend. Use this tool if you want to customize analysis and visualization of lipid-protein interactions and want to scale-up your workflow beyond the capabilities of the <a href="https://prolint.ca" target="_blank">ProLint webserver</a>. 

To get familiar with **prolintpy** please read the <a href="https://prolint.github.io/prolintpy" target="_blank">documentation</a>. 

## Installation 

Getting `prolintpy` is quite simple, especially if you have `conda` installed. 

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

