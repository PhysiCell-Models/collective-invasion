# PhysiCell: Extracellular Matrix Modeling

An extracellular matrix (ECM) model based on leader-follower (K14+/K14-) interactions created using PhysiCell.


## Overview:

The ECM model is built using [PhysiCell](https://github.com/MathCancer/PhysiCell), an agent-based, multicellular 3-D modeling framework written in C++.  If you are not already familiar with PhysiCell then you should do so before looking at the code or trying to modify it.  The ECM model has three components: `density` (thickness of the fibers), `anisotropy` (agreement of alignment) and `fiber alignment` (three coordinate vector of overall fiber direction).  There are two different cell types in the model: leader cells (K14+) and follower cells (K14-).  In the model, leader cells can influence the ECM but are not affected by the ECM while follower cells are effected by the ECM but cannot change it.

The philosophy behind the model is that ECM `density` (from 0-1) remains constant, the `fiber alignment` vector gives the overall alignment of the fiber and `anisotropy` (0-1) represents how aligned the ECM fibers are.  An `anisotropy` value closer to 0 mean the fibers are less aligned and values closer to 1 mean the ECM fibers are more aligned.  The `fiber alignment` vector represents the average fiber direction.  As leader cells move through the ECM, they align the fibers in the direction of their moment, increasing the anisotropy in proportion to the cells speed.  Follower cells "read" the ECM and use the fibers to influence their direction.  Density influences cell speed as the cell gets stuck if the fibers are too thick and if the fibers are too thin then the cells have nothing to "grab".  The fiber direction influences the motility vector of the cells and the amount of influence is controlled by the anisotropy.

The model was developed by Paul Macklin, John Metzcar, Brandon Fischer, Daniel Murphy and Ben Duggan.  Significant help was provided by Randy Heiland and xxx at Johns Hopkins University.


**Reference:** A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: [10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)

**Latest info:**  follow [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)


## Install

Before trying to run the code, you should ensure that your computer can run PhysiCell.  You should read the PhysiCell [quick start guide](https://github.com/MathCancer/PhysiCell/blob/master/Quickstart.pdf) and ensure that your computer has the proper dependence installed.

It is recommended that you download the model as [release](https://github.com/MathCancer/AMIGOS-invasion/releases) of the code.  However, you could download or clone the master branch for the most up to date, stable version of the code.


## How to run

This model runs similar to most other PhysiCell projects, however, there are some differences that are noted below.  

After navigate to the root directory, run `make`.  This will compile the PhysiCell code and create an executable names `AMIGOS-invasion`.  If you have issues compiling or running the code you should consult the PhysiCell [quick start guide](https://github.com/MathCancer/PhysiCell/blob/master/Quickstart.pdf).

### Key makefile rules:

```
make                    : compiles the ECM Modeling code
make clean              : removes all .o files and the executable, so that the next `make` recompiles the entire project 
make data-cleanup       : clears out all simulation data
make data-cleanup-all   : clears out all simulation data (*.xml, *.svg, *.pov, Output/*, SVG/*)
make data-cleanup-light : lighly cleans out the simulation data (*.xml, *.mat, *.svg *.pov, Output/*.mat, Output/*.xml, Output/*.svg, Output/*.png, SCG/*)
```

### Settings

Most of the model settings can be edited in the [PhysiCell_settings.xml](config/PhysiCell_settings.xml).  While you don't have to edit this file you can get different behavior by changing the settings in the file.  Below are some of the settings you may want to play with.

### Outputs

The outputs from the simulation can be found in the [Output](Output/) and [SVG](SVG/) folders.  Leader (K14+) cells are drawn as blue circles and follower (K14-) cells are drawn as yellow circles.

In addition to the standard PhysiCell outputs, our model outputs an ECM specific MATLAB file at each save time.  The file is saved to the [Output](Output/) folder with the form `outputxxxxxxxx_ECM.mat` where `xxxxxxxx` is the current time step.  The MATLAB file contains one array named `ECM_Data` with each column representing a voxel (indexed using voxel id) and rows representing x voxel cord coordinate, y voxel coordinate, z voxel coordinate, ECM anisotropy, ECM density, ECM x alignment, ECM y alignment, ECM z alignment, oxygen gradient x, oxygen gradient y, oxygen gradient z from index 1 to 11.  

This data can be visualized using the scripts in [python_imaing](python_imaging/), particularly [finished-combined-plot.py](python_imaging/finished-combined-plot.py).  This script will make a plot with cells (leaders are blue and followers are yellow), a contour plot showing the oxygen gradient and a quiver plot showing the fiber alignment, scaled by anisotropy.


## Release summary:

* ECM model for UROC poster

See [changes.md](changes.md) for the full change log. 
 
