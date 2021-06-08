<h3 align="center"><img src="https://i.imgur.com/HomQmrp.png" alt="logo" height="200px"></h3>
<p align="center">A python package for the automated analysis and visualization of lipid-protein interactions.</p>


<p align="center">
<a href="https://mybinder.org/v2/gh/ProLint/tutorials/main"><img src="https://mybinder.org/badge_logo.svg"></a>
<a href="https://gitter.im/ProLint/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge"><img src="https://badges.gitter.im/ProLint/community.svg"><a/>
<a href="./LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
<a href="https://github.com/ProLint/prolintpy/releases"><img src="https://img.shields.io/github/v/release/ProLint/prolintpy.svg"></a>
  
<!--   [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ProLint/tutorials/main) -->
</p>


prolintpy is a lightweight python library that is used by the ProLint webserver on the backend to analyze **Pro**tein-**L**ipid **int**eractions.. Use this tool if you want to customize analysis and visualization of lipid-protein interactions and want to scale-up your workflow beyond the capabilities of the <a href="https://prolint.ca" target="_blank">ProLint webserver</a>. 

To get familiar with **prolintpy** please read the <a href="https://prolint.github.io/prolintpy" target="_blank">Documentation</a>. You can also launch the binders: 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ProLint/tutorials/main) and test **prolintpy** without having to install it. 
We provide several tutorials in the form of these binders that you can complete in your browser. Please note, however, that these tutorials are still work in progres and we will improve them during the coming days. If you would like to download the tutorials they are located here: https://github.com/ProLint/tutorials

## Installation 

To install `prolintpy` simply execute: 

```sh
python -m pip install prolintpy
```

This should work on most systems.
If you encounter problems during installation, try the more verbose method below: 

```sh
# create new environment
conda create --name prolint python=3.7
conda activate prolint
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

To import `prolintpy` in your code: 

```python 
# To improve performance of visualization apps, we load everything at once 
# so initial import may take a few seconds to finish. 
import prolintpy as pl
```
Please follow the instructions provided in the <a href="https://prolint.github.io/prolintpy" target="_blank">documentation</a> to get started. Note that, to use the visualization interface of prolintpy, 
you should use JupyterLab. At the top of your notebook file, make sure to call the `output_notebook` function: 

```python
from bokeh.io import output_notebook
output_notebook()
```

Additionally, if you want to use the `show_contact_projection` function, make sure that your installation of `nglview` is working properly. 
Follow the instruction provided <a href="https://github.com/nglviewer/nglview" target="_blank">there</a> to ensure your installation is running correclty.

## Input file requirements
Before you load the data to prolintpy make sure to first remove water & ions beads. Leave only membrane and protein beads in the system. 

## Citation
prolintpy is research software. If you make use of it in work which you publish, please cite it. The BibTeX reference is

```
@article{10.1093/nar/gkab409,
    author = {Sejdiu, Besian I and Tieleman, D Peter},
    title = "{ProLint: a web-based framework for the automated data analysis and visualization of lipid–protein interactions}",
    journal = {Nucleic Acids Research},
    year = {2021},
    month = {05},
    issn = {0305-1048},
    doi = {10.1093/nar/gkab409},
    url = {https://doi.org/10.1093/nar/gkab409},
    note = {gkab409},
    eprint = {https://academic.oup.com/nar/advance-article-pdf/doi/10.1093/nar/gkab409/38270894/gkab409.pdf},
}
```
  
