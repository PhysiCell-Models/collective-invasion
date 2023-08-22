# Collective invasion and extracellular matrix modeling in an cell-based modeling framework

### A model of collective invasion based on multicellular communication mediated through a novel extracellular matrix (ECM) model implemented in PhysiCell.

## Overview
<p align="center"><img alt="NA" src="https://github.com/PhysiCell-Models/collective-invasion/raw/master/Figure-5_Collective_migration_is_robust.png" style="width: 400px; height: 300px;" /></a></p>

The ECM model has three components: `anisotropy` (agreement of alignment), `fiber orientation` (overall orientation of the fibers), and `density` (relative volume filling of ECM fibers).  The collective invasion models uses two different cell types: leader cells and follower cells. Leader cells move up chemotactic gradients and signal their paths by remodeling the ECM. Follower cells alter their motility in response to signals in the ECM as well as chemotaxing when on remodeled ECM and move randomly otherwise. The coupling of ECM signal generating phenotype with ECM signal reading phenotype enables collective behavior (both stigmery and collective invasion) in a range of parameter values. The model and results are explained in depth in [2](#references). See [3](#references) and [4](#references) for a partial biological background. 

The cell and ECM models are built using [PhysiCell](https://github.com/MathCancer/PhysiCell) v 1.4.1, an open source cell-based, multicellular 3-D modeling framework written in C++.  If you are not already familiar with PhysiCell please start by reviewing it (see [1](#references) and [Install](#install) section below) prior to diving deeply into the modeling code here.  

A cloud-based, executable version of this model is available at [nanoHUB](https://nanohub.org/tools/physicellecm). A free login is required to access this resource. 

### ECM and leader-follower models and emergent results

The ECM is conceived of as being composed of a set of small ECM units, the properties of which can be described by 3 components representing the average value of each component over the small unit. These ECM components impact the motility of ECM sensitive cells and in turn can be altered by cells. 

- `Density` (a scalar ranging from 0-1) represents how much overall fiber is present relative to how much space is available to be filled with fibers
- `Fiber orientation` (numerically given as a unit vector) represents the overall (mean) fiber orientations
- `Anisotropy` (0-1) represents how aligned the ECM fibers. 

An `anisotropy` value closer to 0 means the fibers are less aligned, producing little signal, and values closer to 1 indicatd highly aligned fibers that produce a strong signal. `Density` influences cell speed as the cell gets stuck if the fibers are too thick and if the fibers are too sparse then the cells have nothing to "grab".  The `fiber orientation` influences the motility vector of the cells and the amount of influence is controlled by the `anisotropy`. ECM modifying cells can align the fibers in the direction of their movement (changing `fiber orientation`), increase the `anisotropy` in proportion to the cell's speed, and modify `density` up or down to a target density. 

In our example, leader cells are endowed with ECM modification abilities while follower cells "read" the ECM with accompanying changes in cell motility. Instantiating both cell phenotypes in one simulation and altering the rates of ECM modification and ratios of cell speed to cell-cell adhesion produce a range of multicellular behaviors - including stigmergy, collective invasion, uncoupled behavior of the two populations, or a homestatis like pattern. 

The models and computational experiments were conceived of and developed by John Metzcar, Ben Duggan, Brandon Fischer, Matthew Murphy and Paul Macklin.  Significant help was provided by Randy Heiland and the overall AMIGOS Team. Work was funded by a joint (AMIGOS) JKTGF and BCRF grant.

## Installing PhysiCell

PhysiCell has a small number of dependencies. Before trying to run the code, you should ensure that your computer can run PhysiCell. Read the PhysiCell [Quick Start Guide](https://github.com/MathCancer/PhysiCell/blob/master/documentation/Quickstart.md) to begin working with PhysiCell and ensure that your system has the proper dependencies installed or follow the [Mac install](https://www.youtube.com/watch?v=Sq9nfKS5U0E&list=PL1fyIV-yPAYzzOVxfGsL90a5KTTh8gSW2&index=2) or [Windows install](https://www.youtube.com/watch?v=hIP4JUrViRA) videos. Note that the videos include more than is required to run the ECM and collective invasion models. Following the install videos through around the half way point of each video is sufficient. SBML and the PhysiCell Model Builder are not required. Python, Matplotlib and other basic Python packages are required to produce visualizations behind the SVGs produced by core PhysiCell. 

Additional install support (and otherwise) can be found be generating an issue at [SourceForge](https://sourceforge.net/projects/physicell/) or on the [PhysiCell Slack workspace](https://join.slack.com/t/physicellcomm-sf93727/shared_invite/zt-qj1av6yd-yVeer8VkQaNDjDz7fF00jA). 

Note that we recommend downloading the model [release](https://github.com/MathCancer/AMIGOS-invasion/releases). However, you could download or clone the main branch for the most up to date, stable version of the code.

## How to run

This model runs similarly to most other PhysiCell projects, however, there are some differences that are noted below.  

After navigating to the root directory, run `make`.  This will compile the PhysiCell code and create an executable named `AMIGOS-invasion`.  If you have issues compiling or running the code, begin by consulting the PhysiCell [Quick Start Guide](https://github.com/MathCancer/PhysiCell/blob/master/documentation/Quickstart.md) and other documentation.

Please submit issues specific to the collective invasion and ECM models as an [issue](https://github.com/PhysiCell-Models/collective-invasion/issues) on this repo.

## Settings and exploring the collective invasion model

Many model settings can be edited in the [PhysiCell_settings.xml](config/PhysiCell_settings.xml).  While you don't have to edit this file you can get different behavior by changing the settings in the file.  

You can follow these suggestions to familiarize yourself with the model.

1) Run this model using the default parameters. This will be produce the right hand side of the figure at the beginning of this README (Figure 5 of Reference [2]). 

2) Change anisotropy_increase_rate to 0.001 and fiber_realignment_rate to 1. This will decrease the collective behavior - eliminating collective invasion and leaving many followers in the center of the domain.

3) Change discrete_ECM_remodeling to 0. Note now the ECM is instantly remodeled - generating strong, clear signals for followers to read (view the anisotropy field). This will recover behavior similar to the default parameters and produce the left side of the figure included above.

4) Change default_cell_speed to 1.0 and then 0.25 - producing stigmergy and a non-changing morphology.

5) Now, you can enjoy yourself changing other parameters and creating new responses.

Note that code does not need recompiled in between parameter changes; the executable will parse the changes to the xml. 

## ECM Specific Outputs

In addition to the standard PhysiCell outputs, our model outputs an ECM specific MATLAB file at each save time.  The file is saved to the `output` folder with the form `outputxxxxxxxx_ECM.mat` where `xxxxxxxx` is the _i_ th simulation output.  The MATLAB file contains one array named `ECM_Data` with each column representing a voxel (indexed using voxel id) and rows representing the x, y, and z voxel coordinates, ECM anisotropy, ECM density, and fiber orientation x, y, and z components.  

This data can be visualized using the scripts in [python_imaging](python_imaging/). We provide general image production through a general template script [image_processing_script.py](python_imaging/image_processing_script.py) which accesses the *PhysiCellPlotter* class in the module *Image processing for PhysiCell* in [image_processing_for_physicell.py](python_imaging/image_processing_for_physicell.py). For image production settings optimzied for the default parameter settings, see [partial_history_multilevel_contour_still.py](python_imaging/partial_history_multilevel_contour_still.py) and [partial_history_multilevel_contour_movie.py](python_imaging/partial_history_multilevel_contour_movie.py). This script will make an overlaying composite plot of cells (leaders are blue and followers are yellow), a contour plot showing oxygen (in red) and a quiver plots showing cell movement history.

## Running PhysiCell simulations across a team

PhysiCell can be run in a distributed fashion across a team using DAPT: Distributed Automated Parameter Testing [5]. See the coder [here](https://github.com/BenSDuggan/DAPT) including a [detailed PhysiCell example](https://github.com/PhysiCell-Tools/DAPT-example).

## Some key makefile rules

```
make                    : compiles the simple test codes
make                    : compiles the leader-follower model
make fibrosis           : compiles the fibrosis model
make invasive_carcinoma : compiles the invasive carcinoma model
make clean              : removes all .o files and the executable, so that the next `make` recompiles the entire project 
make data-cleanup       : clears out all simulation data
make data-cleanup-all   : clears out all simulation data (*.xml, *.svg, *.pov, Output/*, SVG/*)
make data-cleanup-light : lightly cleans out the simulation data (*.xml, *.mat, *.svg *.pov, Output/*.mat, Output/*.xml, Output/*.svg, Output/*.png, SCG/*)
make zip-source         : compresses all files required to reproduce simulation as well as the python imaging folder
```

See makefile for additional rules. 

## References

[1] Ghaffarizadeh, A, Heiland, R, Friedman, SH, Mumenthaler, SM, and Macklin, P. PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: [10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)

[2] Metzcar, J, Duggan, BS, Fischer, B, Murphy, M, Macklin, P. A novel model of multicellular communication through extracellular matrix microstructure. bioRxiv 2022.11.21.514608; doi: [https://doi.org/10.1101/2022.11.21.514608](https://doi.org/10.1101/2022.11.21.514608)

[3] Cheung KJ, Gabrielson E, Werb Z, Ewald AJ. Collective invasion in breast cancer requires a conserved basal epithelial program. Cell 2013; 155(7):1639-51. [10.1016/j.cell.2013.11.029](10.1016/j.cell.2013.11.029)

[4] Nguyen-Ngoc KV\*, Cheung KJ*, Brenot A, Shamir ER, Gray RS, Hines WC, Yaswen P, Werb Z, Ewald AJ. The ECM microenvironment regulates collective migration and local dissemination in normal and malignant mammary epithelium. Proceedings of the National Academy of Science 2012; [10.1073/pnas.1212834109](10.1073/pnas.1212834109) *Co-First Authors. PMCID: PMC3465416

[5] Duggan, BS, Metzcar, J, and Macklin, P (2021). DAPT: A package enabling distributed automated parameter testing. Gigabyte 2021, 1–10. [10.46471/gigabyte.22](10.46471/gigabyte.22).


**Latest PhysiCell info:**  follow [@PhysiCell](https://twitter.com/PhysiCell) on Twitter (http://twitter.com/PhysiCell)
 
