# Extracellular matrix modeling in an cell-based modeling framework

An extracellular matrix (ECM) model implemented in PhysiCell.


## Overview:

The ECM model is built using [PhysiCell](https://github.com/MathCancer/PhysiCell), an open source cell-based, multicellular 3-D modeling framework written in C++.  If you are not already familiar with PhysiCell please start by reviewing it (see references and **Install** section below) prior to diving deeply into the modeling code here.  

The ECM model has three components: `density` (relative quantity of the fibers), `anisotropy` (agreement of alignment) and `fiber direction` (vector of overall fiber direction).  In the example code, we use two different cell types in the model: leader cells (K14+) and follower cells (K14-). Leader cells can influence the ECM but are not affected by the ECM while follower cells are affected by the ECM but cannot change it.

The ECM is conceived of as being composed of a set of small ECM units, the properties of which can be described by 3 components representing the average value of each component over the small unit. These ECM components are capable of impacting the motility of ECM sensitive cells and in turn can be altered by cells. The `density` (a scalar ranging from 0-1) represents how much overall fiber is present relative to how much space is available to be filled with fibers, the `fiber direction` (a unit vector) represents the overall (mean) fiber direction, and `anisotropy` (0-1) represents how aligned the ECM fibers are as well as elements of fiber resistance to deformation and realignment.  An `anisotropy` value closer to 0 means the fibers are less aligned and values closer to 1 mean the ECM fibers are more aligned and more resistant to realignment. `Density` influences cell speed as the cell gets stuck if the fibers are too thick and if the fibers are too sparse then the cells have nothing to \`\`grab".  The `fiber direction` influences the motility vector of the cells and the amount of influence is controlled by the `anisotropy`. ECM modifying cells can align the fibers in the direction of their movement (changing `fiber direction`), increase the `anisotropy` in proportion to the cell's speed, and modify `density` up or down to a target density. In our example, leader cells are endowed with ECM modification abilities while follower cells \`\`read" the ECM with accompanying changes in cell motility. 

The model was conceived of and developed by John Metzcar, Ben Duggan, Brandon Fischer, Daniel Murphy and Paul Macklin.  Significant help was provided by Randy Heiland and the overall AMIGOS Team. Work was funded by a joint (AMIGOS) JKTGF and BCRF grant.


**References:** 

_Main code references_: 

Ghaffarizadeh, A, Heiland, R, Friedman, SH, Mumenthaler, SM, and Macklin, P. PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: [10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)

Metzcar, J, Duggan, B, Fisher, B, Murphy, D, Macklin, P. A simple extracellular matrix model for an agent based modeling framework, 2019. Pending submission.

_Biology references_:

Cheung KJ, Gabrielson E, Werb Z, Ewald AJ. Collective invasion in breast cancer requires a conserved basal epithelial program. Cell 2013; 155(7):1639-51.

Nguyen-Ngoc KV\*, Cheung KJ*, Brenot A, Shamir ER, Gray RS, Hines WC, Yaswen P, Werb Z, Ewald AJ. The ECM microenvironment regulates collective migration and local dissemination in normal and malignant mammary epithelium. Proceedings of the National Academy of Science 2012; 10.1073/pnas.1212834109 *Co-First Authors. PMCID: PMC3465416




**Latest info:**  follow [@PhysiCell](https://twitter.com/PhysiCell) on Twitter (http://twitter.com/PhysiCell)


## Install

Before trying to run the code, you should ensure that your computer can run PhysiCell.  Read the PhysiCell [Quick Start Guide](https://github.com/MathCancer/PhysiCell/blob/master/Quickstart.pdf) to begin working with PhysiCell and ensure that your system has the proper dependencies installed.

It is recommended that you download the model as a [release](https://github.com/MathCancer/AMIGOS-invasion/releases) of the code.  However, you could download or clone the master branch for the most up to date, stable version of the code.


## How to run

This model runs similarly to most other PhysiCell projects, however, there are some differences that are noted below.  

After navigating to the root directory, run `make`.  This will compile the PhysiCell code and create an executable named `AMIGOS-invasion`.  If you have issues compiling or running the code, begin by consulting the PhysiCell [Quick Start Guide](https://github.com/MathCancer/PhysiCell/blob/master/Quickstart.pdf) and other documentation.

### Key makefile rules:

```
make                    : compiles the ECM Modeling code
make clean              : removes all .o files and the executable, so that the next `make` recompiles the entire project 
make data-cleanup       : clears out all simulation data
make data-cleanup-all   : clears out all simulation data (*.xml, *.svg, *.pov, Output/*, SVG/*)
make data-cleanup-light : lightly cleans out the simulation data (*.xml, *.mat, *.svg *.pov, Output/*.mat, Output/*.xml, Output/*.svg, Output/*.png, SCG/*)
```

### Settings

Most of the model settings can be edited in the [PhysiCell_settings.xml](config/PhysiCell_settings.xml).  While you don't have to edit this file you can get different behavior by changing the settings in the file.  Below are some of the settings you may want to play with.

Leader and follower adhesion levels + default cell speeds - values of 6.5 and 10 for leader and follower adhesion (letting leader adhesion = follower adhesion) combined with default cell speeds of 0.3 and 0.5 respectively produce interesting simulations displaying leader-follower behavior. 

### Outputs

The outputs from the simulation can be found in the [Output](Output/) and [SVG](SVG/) folders.  Leader (K14+) cells are drawn as blue circles and follower (K14-) cells are drawn as yellow circles.

In addition to the standard PhysiCell outputs, our model outputs an ECM specific MATLAB file at each save time.  The file is saved to the [Output](Output/) folder with the form `outputxxxxxxxx_ECM.mat` where `xxxxxxxx` is the current time step.  The MATLAB file contains one array named `ECM_Data` with each column representing a voxel (indexed using voxel id) and rows representing the x, y, and z voxel coordinates, ECM anisotropy, ECM density, fiber direction x, y, and z components, and finally the oxygen gradient x, y, and z components from index 1 to 11 (assuming the file is being read into Matlab).  

This data can be visualized using the scripts in [python_imaging](python_imaging/), particularly [finished-combined-plot.py](python_imaging/finished-combined-plot.py).  This script will make an overlaying composite plot of cells (leaders are blue and followers are yellow), a contour plot showing the oxygen gradient and a quiver plot showing the fiber alignment, scaled by anisotropy.


## Release summary:

* ECM model for UROC poster

[//]:# (See [changes.md](changes.md) for the full change log.)
 
