# About

[ProLint](https://www.prolint.ca) is a framework dedicated to the analysis and visualization of lipid-protein interactions. <br>
ProLint is not just a single tool or collection of scripts, but a list of dedicated programs and software libraries.

Briefly, it includes:
<ul>
<li>The ProLint webserver</li>
<li>The ProLint database</li>
<li>The g_surf program</li>
<li>The stand-alone python package prolintpy</li>
</ul>

The following documentation relates to the stand-alone python package `prolintpy`.

ProLint is also a webservice in that it can store results and allow users world-wide
to interact with the data regardless of their knowledge of Molecular Dynamics (MD) Simulations or coding experience.
Hosting and sharing of data files is done free of charge.

The software that is made available under ProLint is either a binary (that does not require installation) or comes with a very
easy to install process. All of the code is shared and made available on [GitHub](https://github.com/ProLint/ProLint) with open-source licenses. <br>
Interested users can download the code and use all of the software and programs locally.

## Goals and Objectives
ProLint as a webserver and webservice has four core objectives. These are listed below and explained in more detail in this section.

<ol>
<li>Modular Analysis.</li>
<li>Interactive Visualization.</li>
<li>Automation and Scalability.</li>
<li>Accessibility and Shareability.</li>
</ol>

### Modular Analysis

Currently, the time it takes to analyze and interpret data usually exceeds the time
to set up and perform Molecular Dynamics (MD) Simulations. Add to this the
continuing trend of increased system sizes, simulated for longer time scales and
the need for a much bigger set of proteins, it becomes clear how analysis and
visualization of protein-lipid contacts can become a significant bottleneck.
Analysis is always context-dependent and users will have varying levels of interest for different aspects of the same data.
ProLint recognizes this and implements several analysis protocols to account for it.


ProLint
bridges this gap between data generation and obtaining insight about relevant interactions.
To learn how ProLint analyses contacts, see the [Analysis](analysis.md) and [Metrics](metrics.md) sections.

### Interactive Visualization

ProLint comes with a built-in visualization library that quickly and effortlessly displays
protein-lipid interactions. It uses modern visualization libraries like bokeh, d3.js,
nglviewer to do the heavy-lifting and allows the user to quickly see the presence (or absence) of
protein-lipid contacts in the simulated system.
The line between analysis and visualization becomes indistinguishable once users have the ability to
interact with the data being visualized.
The design of all visualization applications by ProLint is centered around user interactivity!
To learn how ProLint visualizes data, see the [Visualization](visualization.md) sections.

### Automation and Scalability

It will become increasingly important to automate the generation of protein-lipid contacts for
any arbitrary protein. Such databases will be of immense importance to the scientific community.
The library presented here `prolintpy` is perfectly suited to be used as part of pipelines to automate tasks and scale-up workflow.


## Accessibility and Shareability

Papers studying protein-lipid contacts are closed systems. They use widely different analyses
protocols and visualization tools that make it difficult to compare results across different
studies. Getting the original data that was used in published works is also a challenge.
ProLint allows for contact information to be extracted from simulated systems and stored
in small files that can be easily shared and made publicly available.
Web-based visualization and cloud-based hosting allow ProLint to share MD results with the whole scientific community.
Thanks to its intuitive visualization applications that require no coding knowledge or MD experience, the benefit of the shared data extends far beyond the computational community.


