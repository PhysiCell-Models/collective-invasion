# A simple framework for agent-based modeling with extracellular matrix

<!-- ### A simple framework for agent-based modeling with extracellular matrix -->

## Framework overview
<p align="center"><img alt="NA" src="https://github.com/PhysiCell-Models/collective-invasion/raw/master/ECM_framework_schematic.png" style="width: 400px; height: 300px;" /></a></p>

In this framework [[1](#references)] for modeling the ECM and cell-ECM interactions, we divide the ECM into volumetric elements that track local ECM density, alignment, and overall anisotropy (local microstructure). Individual cell agents can locally remodel each of these properties, while these properties can in turn influence cell behavior including changes in migration speed, chemotactic response, ECM contact guidance, proliferation, death, secretion, and differentiation. It is implemented as an extension of the open source package PhysiCell [[2]](#references) and can readily be used to incorporate local ECM effects into agent-based models, in particular through PhysiCell rules [[3]](#references). 

For additional details and background on the components of this framework, see [[1](#references)]

## Overview of repository structure and key files

### ECM and cell-ECM core files

All the examples use the same core source code files. These are the core of the framework and its extension to PhysiCell. They are:

- `extracellular_matrix.h/extracellular_matrix.cpp`
    - Contains the ECM element and ECM mesh class definitions and other initilization routine
- `cell_ECM_interactions.h/cell_ECM_matrix.cpp`
    - Contains functions for the bidirectional cell-ECM interactions, generating a default ECM compatible agent, and custom output routines

These files are currently located in `custom_modules`. 

### Model files

There are three example models, several variants of the leader-follower model, and simple test models. See the below `Compling and running sample models` below for details on those models. 

There is one makefile for all compilation. It is located in the `root`. 

Each base model has a main_XXX.cpp (in `root`), custom source code (in `custom_code`), and one to several additional model specficifcation files in `config` - always including an XML-based model config file and sometimes including PhysiCell rules and initial cell position files (both .csv's). The config directory also contains special XML-bassed model files for testing and making stochastic replicates. 


## Using the framework

The cell-ECM interaction framework is built as an extension to [PhysiCell](https://github.com/MathCancer/PhysiCell) v 1.12, an open source cell-based, multicellular 3-D modeling framework written in C++. If you are not already familiar with PhysiCell we suggest begining by reviewing it (see [2](#references)). Note, that the ECM framework stands alone - coming with all necessary source code to compile and execute the sample models and develop new models. However, as the ECM framework is a direct extension of PhysiCell, we suggest using the PhysiCell install guides for your system - with current guides [here](https://github.com/physicell-training/ws2023/blob/main/agenda.md) and more generally at the [PhysiCell training repository](https://github.com/physicell-training). 

### Compiling and running sample models

There are 3 main sample models as well as a series of simple tests. Listed below are instructions for making and running each model and its variants, grouped by which executable is required. 

#### _Simple tests_:

`make` - compiles the AMIGOS-invasion executable

The following output directories will need made to run the simulations below: `simple_test0`, `simple_test1`, `simple_test3`, and `simple_test4`. `simple_test2` is already in the repository.

Simple tests - demonstrating the main cell-ECM interactions as one way (either ECM remodeling or ECM following) experiments:

- `./AMIGOS-invasion config/simple_test0_straight_ECM.xml` - ECM following: random 1-D motion along vertically oriented ECM
- `./AMIGOS-invasion config/simple_test1_cell_march.xml` - ECM remodeling: realignment of randomly oriented fiber orientations (Figure 2a from [[1](#references)])
- `./AMIGOS-invasion config/simple_test2_random_1_D_circles.xml` - ECM following: random 1-D motion along circularly oriented ECM (Figure 2b from [[1](#references)])
- `./AMIGOS-invasion config/simple_test3_directed_circular_motion.xml` - ECM following: combining a second direction, a chemotactic gradient, with ECM following on circularly oriented ECM (Figure 2c from [[1](#references)])
- `./AMIGOS-invasion config/simple_test4_split_ECM.xml` - ECM following: combining a second direction, a chemotactic gradient, with ECM following on split, diagnoally oriented ECM (SM Figure 3 from [[1](#references)])

#### _Wound headling and fibrosis model_:

`make fibrosis` - compiles fibrosis executable

- `./fibrosis config/fibrosis.xml` - simulated tissue insult is cleared by macrophages, which recruit fibroblasts that increase ECM density in the presence of macrophages, leading to hyperdense ECM that is relatively impenetrable cells surrounding a region of relatively less dense ECM where the tissue insult occurred (Figure 3 from [[1](#references)]).

Output will go to `fibrosis_test`

#### _Basement membrane degradation and stromal invasion_:

`make invasive_carcinoma` - compiles invasive carcinoma executable

- `./invasive_carcinoma config/invasive_carcinoma.xml` - simulation of basement membrane degradation by tumor recruited fibroblasts, leading to invasion of stroma by previously _in situ_ tumor (Figure 4 from [[1](#references)]).

Output will go to `invasive_carcinoma_output`

#### _Leader-follower and collective migration model_:

`make` - compiles the AMIGOS-invasion executable

The following output directories will need made to run the simulations below: `adh_0_repulsion_0_speed_10_no_reading`, `adh_0_repulsion_0_speed_10_no_writing`, `adh_0_repulsion_0_speed_10`, `adh_10_replusion_25_speed_080`, `adh_10_replusion_25_speed_050`, `adh_10_replusion_25_speed_010`, and `leader_follower_decreased_remodeling`. `leader_follower` is already in the repository. 

- `./AMIGOS-invasion writing_only.xml`- No contact guidance in the follower cell population - notably lacks outward migration of fiber following cells ("followers") (Figure 5a from [[1](#references)])
- `./AMIGOS-invasion reading_only.xml` - No production of directional cues in ECM for followers to follow - notably lacks outward migration of fiber following cells (Figure 5b from [[1](#references)])
- `./AMIGOS-invasion writing_and_reading.xml` - Enables the two cell populations (leader and follower), with leaders producing paths in the ECM that followers follow (stigmergy) (Figure 5c from [[1](#references)])
- `./AMIGOS-invasion leader_follower_instant_speed_080.xml` - Produces stigmery (Figure 6a from [[1](#references)])
- `./AMIGOS-invasion leader_follower_instant_speed_050.xml` - Produces leader-follower collective migration (Figure 6b from [[1](#references)])
- `./AMIGOS-invasion leader_follower_instant_speed_010.xml` - Produces no change in cell arrangement/pattern (Figure 6c from [[1](#references)])
- `./AMIGOS-invasion leader_follower_model.xml` - Demonstrates leader-follower collective migration, even in the case of non-instant ECM remodeling (Figure 7b from [[1](#references)])
- `./AMIGOS-invasion leader_follower_model_decreased_remodeling_rates.xml` - Demonstrates that leader-follower collective migration is relatively senstive to remodeling rate parameters (SM Figure 4b from [[1](#references)])

### Running through the Studio

All the above models (and in general any PhysiCell model) can be run via the PhysiCell Studio ([User Guide](https://github.com/PhysiCell-Tools/Studio-Guide/blob/main/README.md) and reference [[4]](#references)). 

Briefly, to run an ECM-based model AND enable visualization of the ECM variables (anisotropy, density, and orientation), use the following pattern to start the Studio:

`python path_to_studio_directory/studio_ecm.py -e [executable_name] -c [config/config_file_name]`

This assumes you are invoking python in the same directory as the executable and that your config file is in `config`. Note that in the command above, we use `studio_ecm.py` NOT `studio.py`. Using `studio_ecm.py` will add built-in ECM field visualization to the Studio. However, this is not an officially supported feature of the Studio, so it may at some point be deprecated without notice. We will attempt to maintain compatiability. 

Note that the ECM-based models are modifiable in the regular version of the studio, but the  ECM visualization is not. Please see the PhysiCell Studio Guide and preprint for general information and details on the Studio.

### nanoHUB

A cloud-based, executable version of the leader-follower model is available at [nanoHUB](https://nanohub.org/tools/physicellecm). A free login is required to access this resource. 


## ECM Specific Outputs

In addition to the standard PhysiCell outputs, our model outputs an ECM specific MATLAB file at each save time.  The file is saved to the `output` folder with the form `outputxxxxxxxx_ECM.mat` where `xxxxxxxx` is the _i_ th simulation output.  The MATLAB file contains one array named `ECM_Data` with each column representing a voxel (indexed using voxel id) and rows representing the x, y, and z voxel coordinates, ECM anisotropy, ECM density, and fiber orientation x, y, and z components.  

ECM anisotropy and density can also be output to non-diffusing fields through the function `copy_ECM_data_to_BioFVM` in `cell_ECM_interactions.cpp`. Currently, this is enabled by uncommenting the line `// copy_ECM_data_to_BioFVM();` in each models main file. Note that to use this feature, the model needs to have the fields `ECM_anisotropy` and `ECM_density` in the model config file (xml file) and the ECM element size has to match the diffusion voxel size. These fields can be added to any model config file through the Studio - make a new field, change the field name to match the above exactly and accept the defaults (0 for everything). Do this twice. 

## ECM visualization

This data can be visualized using the README and scripts in [python_imaging](python_imaging/). We provide general image production through a general template script [image_processing_script.py](python_imaging/image_processing_script.py) which accesses the *PhysiCellPlotter* class in the module *Image processing for PhysiCell* in [image_processing_for_physicell.py](python_imaging/image_processing_for_physicell.py). For image production settings optimzied for the default parameter settings, see [partial_history_multilevel_contour_still.py](python_imaging/partial_history_multilevel_contour_still.py) and [partial_history_multilevel_contour_movie.py](python_imaging/partial_history_multilevel_contour_movie.py). This script will make an overlaying composite plot of cells (leaders are blue and followers are yellow), a contour plot showing oxygen (in red) and a quiver plots showing cell movement history. Note that we do our best to ensure that all code and scripts in `python_imaging` work without alteration and as expected - and please consider it to be a preliminary release that is not guanteed to work and that may change in the future. 

For rapid visualization, the Studio can be used - in either regular or ECM mode. Note that if ECM anisotropy and density are output to the standard microenvironment outputs, they can visualized with regular version of the Studio. 

## Running PhysiCell simulations across a team

PhysiCell can be run in a distributed fashion across a team using DAPT: Distributed Automated Parameter Testing [5]. See the code [here](https://github.com/BenSDuggan/DAPT) including a [detailed PhysiCell example](https://github.com/PhysiCell-Tools/DAPT-example).

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

## Future work

- Remove deprecated user_parameters from all model files and code base. 
- Add "Exploration of leader-follower collective migration model" into README (using previous, but currently out of date material). 
    - Could include several model walk throughs by video
    - Could include more in depth context for the leader-follower collective migration
    - Could include more on the ECM model details (again out of date material is available for updating)
- Review optimal way to include `copy_ECM_data_to_BioFVM`
    - Will possibly add an XML parameter for this
    - Will review possibly calling at mechanics or phenotype time step (currently being called at diffusion time step)
- Additonal future work is included in [[1](#references)]

## Acknowledgements

This work was funded in part by a joint (AMIGOS) JKTGF and BCRF grant. We thank Margherita Botticelli for many productive conversations on the cell-ECM interaction code.

## References

[1] Metzcar, J., Duggan, B.S., Fischer, B., Murphy, M., Heiland, R., Macklin, P, 2024. A simple framework for agent-based modeling with extracellular matrix. bioRxiv 2022.11.21.514608; https://doi.org/10.1101/2022.11.21.514608

[2] Ghaffarizadeh, A., Heiland, R., Friedman, S.H., Mumenthaler, S.M., and Macklin, P. PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. https://dx.doi.org/10.1371/journal.pcbi.1005991

[3] Johnson, J.A.I., Stein-O’Brien, G.L., Booth, M., Heiland, R., Kurtoglu, F., Bergman, D.R., Bucher, E., Deshpande, A., Forjaz, A., Getz, M., Godet, I., Lyman, M., Metzcar, J., Mitchell, J., Raddatz, A., Rocha, H., Solorzano, J., Sundus, A., Wang, Y., Gilkes, D., Kagohara, L.T., Kiemen, A.L., Thompson, E.D., Wirtz, D., Wu, P.-H., Zaidi, N., Zheng, L., Zimmerman, J.W., Jaffee, E.M., Chang, Y.H., Coussens, L.M., Gray, J.W., Heiser, L.M., Fertig, E.J., Macklin, P.. Digitize your Biology! Modeling multicellular systems through interpretable cell behavior. bioRxiv, 2023. https://doi.org/10.1101/2023.09.17.557982


[4] Heiland, R., Bergman, D., Lyons, B., Cass, J., Rocha, H.L., Ruscone, M., Noël, V., Macklin, P.. PhysiCell Studio: a graphical tool to make agent-based modeling more accessible. bioRxiv, 2023 https://doi.org/10.1101/2023.10.24.563727

[5] Duggan, B.S., Metzcar, J., and Macklin, P. DAPT: A package enabling distributed automated parameter testing. Gigabyte 2021, 1–10, 2021. https://doi.org/10.46471/gigabyte.22



**Latest PhysiCell info:**  follow [@PhysiCell](https://twitter.com/PhysiCell) on Twitter (http://twitter.com/PhysiCell)