# A simple framework for agent-based modeling with extracellular matrix

<!-- ### A simple framework for agent-based modeling with extracellular matrix -->

## Framework overview
<p align="center"><img alt="NA" src="https://github.com/PhysiCell-Models/collective-invasion/raw/master/ECM_framework_schematic.png" style="width: 400px; height: 300px;" /></a></p>

In this framework [[1](#references)] for modeling the ECM and cell-ECM interactions, we divide the ECM into volumetric elements that track local ECM density, alignment, and overall anisotropy (local microstructure). Individual cell agents can locally remodel each of these properties, while these properties can in turn influence cell behavior including changes in migration speed, chemotactic response, ECM contact guidance, proliferation, death, secretion, and differentiation. It is implemented as an extension of the open source package PhysiCell [[2]](#references) and can readily be used to incorporate local ECM effects into agent-based models, in particular through PhysiCell rules [[3]](#references). 

For additional details and background on the components of this framework, see [1](#references)

## Using the framework

The cell-ECM interaction framework is built as an extension to [PhysiCell](https://github.com/MathCancer/PhysiCell) v 1.12, an open source cell-based, multicellular 3-D modeling framework written in C++. If you are not already familiar with PhysiCell we suggest begining by reviewing it (see [2](#references)). Note, that the ECM framework stands alone - coming with all necessary source code to compile and execute the sample models and develop new models. However, as the ECM framework is a direct extension of PhysiCell, we suggest using the PhysiCell install guides for your system - with current guides [here](https://github.com/physicell-training/ws2023/blob/main/agenda.md) and more generally at the [PhysiCell training repository](https://github.com/physicell-training). 

### Compiling and running sample models

There are 3 main sample models as well as a series of simple tests. Listed below are instructions for making and running each model and its variants, grouped by which executable is required. 

#### _Simple tests_:

`make` - compiles the AMIGOS-invasion executable

Simple tests - demonstrating the main cell-ECM interactions as one way (either ECM remodeling or ECM following) experiments:

- `./AMIGOS-invasion config/simple_test0_straight_ECM.xml` - ECM following: random 1-D motion along vertically oriented ECM
- `./AMIGOS-invasion config/simple_test1_cell_march.xml` - ECM remodeling: realignment of randomly oriented fiber orientations (Figure 2a from [[1](#references)])
- `./AMIGOS-invasion config/simple_test2_random_1_D_circles.xml` - ECM following: random 1-D motion along circularly oriented ECM (Figure 2b from [[1](#references)])
- `./AMIGOS-invasion config/simple_test3_directed_circular_motion.xml` - ECM following: combining a second direction, a chemotactic gradient, with ECM following on circularly oriented ECM (Figure 2c from [[1](#references)])
- `./AMIGOS-invasion config/simple_test4_split_ECM.xml` - ECM following: combining a second direction, a chemotactic gradient, with ECM following on split, diagnoally oriented ECM (SM Figure 3 from [[1](#references)])

#### _Wound headling and fibrosis model_:

`make fibrosis` - compiles fibrosis executable

- `./fibrosis config/fibrosis.xml` - simulated tissue insult is cleared by macrophages, which recruit fibroblasts that increase ECM density in the presence of macrophages, leading to hyperdense ECM that is relatively impenetrable cells surrounding a region of relatively less dense ECM where the tissue insult occurred (Figure 3 from [[1](#references)]).

#### _Basement membrane degradation and stromal invasion_:

`make invasive_carcinoma` - compiles invasive carcinoma executable

- `./invasive_carcinoma config/invasive_carcinoma.xml` - simulation of basement membrane degradation by tumor recruited fibroblasts, leading to invasion of stroma by previously _in situ_ tumor (Figure 4 from [[1](#references)]).

#### _Leader-follower and collective migration model_:

`make` - compiles the AMIGOS-invasion executable

- `./AMIGOS-invasion writing_only.xml`- No contact guidance in the follower cell population - notably lacks outward migration of fiber following cells ("followers") (Figure 5a from [[1](#references)])
- `./AMIGOS-invasion reading_only.xml` - No production of directional cues in ECM for followers to follow - notably lacks outward migration of fiber following cells (Figure 5b from [[1](#references)])
- `./AMIGOS-invasion writing_and_reading.xml` - Enables the two cell populations (leader and follower), with leaders producing paths in the ECM that followers follow (stigmergy) (Figure 5c from [[1](#references)])
- `./AMIGOS-invasion leader_follower_instant_speed_080.xml` - Produces stigmery (Figure 6a from [[1](#references)])
- `./AMIGOS-invasion leader_follower_instant_speed_050.xml` - Produces leader-follower collective migration (Figure 6b from [[1](#references)])
- `./AMIGOS-invasion leader_follower_instant_speed_010.xml` - Produces no change in cell arrangement/pattern (Figure 6c from [[1](#references)])
- `./AMIGOS-invasion leader_follower_model.xml` - Demonstrates leader-follower collective migration, even in the case of non-instant ECM remodeling (Figure 7b from [[1](#references)])
- `./AMIGOS-invasion leader_follower_model_decreased_remodeling_rates.xml` - Demonstrates that leader-follower collective migration is relatively senstive to remodeling rate parameters (SM Figure 4b from [[1](#references)])

### Running through the Studio

All the above models (and in general any PhysiCell model) can be run via the PhysiCell Studio ([User Guide](https://github.com/PhysiCell-Tools/Studio-Guide/blob/main/README.md) and reference [[x]](#references)). 

Briefly, to run an ECM-based model AND enable visualization of the ECM variables (anisotropy, density, and orientation), use the following pattern to start the Studio:

`python path_to_studio_directory/studio_ecm.py -e [executable_name] -c [config/config_file_name]`

This assumes you are invoking python in the same directory as the executable and that your config file is in `config`. Note that in the command above, we use `studio_ecm.py` NOT `studio.py`. Using `studio_ecm.py` will add built in ECM field visualization to the Studio. However, this is not an officially supported feature of the Studio, so it may at some point be deprecated without notice. We will attempt to maintain compatiability. 

Note that the ECM-based models are modifiable in the regular version of the studio, but the ECM visualization is not. Please see the PhysiCell Studio Guide and preprint for general information and details on the Studio.

### nanoHUB

A cloud-based, executable version of the leader-follower model is available at [nanoHUB](https://nanohub.org/tools/physicellecm). A free login is required to access this resource. 


## ECM Specific Outputs

In addition to the standard PhysiCell outputs, our model outputs an ECM specific MATLAB file at each save time.  The file is saved to the `output` folder with the form `outputxxxxxxxx_ECM.mat` where `xxxxxxxx` is the _i_ th simulation output.  The MATLAB file contains one array named `ECM_Data` with each column representing a voxel (indexed using voxel id) and rows representing the x, y, and z voxel coordinates, ECM anisotropy, ECM density, and fiber orientation x, y, and z components.  

## ECM visualization

This data can be visualized using the README and scripts in [python_imaging](python_imaging/). We provide general image production through a general template script [image_processing_script.py](python_imaging/image_processing_script.py) which accesses the *PhysiCellPlotter* class in the module *Image processing for PhysiCell* in [image_processing_for_physicell.py](python_imaging/image_processing_for_physicell.py). For image production settings optimzied for the default parameter settings, see [partial_history_multilevel_contour_still.py](python_imaging/partial_history_multilevel_contour_still.py) and [partial_history_multilevel_contour_movie.py](python_imaging/partial_history_multilevel_contour_movie.py). This script will make an overlaying composite plot of cells (leaders are blue and followers are yellow), a contour plot showing oxygen (in red) and a quiver plots showing cell movement history. Note that we do our best to ensure that all code and scripts in `python_imaging` work without alteration and as expected - and please consider it to be a preliminary release that is not guanteed to work and that may change in the future. 


### ECM and leader-follower models and emergent results


The collective invasion models uses two different cell types: leader cells and follower cells. Leader cells move up chemotactic gradients and signal their paths by remodeling the ECM. Follower cells alter their motility in response to signals in the ECM as well as chemotaxing when on remodeled ECM and move randomly otherwise. The coupling of ECM signal generating phenotype with ECM signal reading phenotype enables collective behavior (both stigmery and collective invasion) in a range of parameter values. The model and results are explained in depth in [2](#references). See [3](#references) and [4](#references) for a partial biological background. 

In our example, leader cells are endowed with ECM modification abilities while follower cells "read" the ECM with accompanying changes in cell motility. Instantiating both cell phenotypes in one simulation and altering the rates of ECM modification and ratios of cell speed to cell-cell adhesion produce a range of multicellular behaviors - including stigmergy, collective invasion, uncoupled behavior of the two populations, or a homestatis like pattern. 



## Settings and exploring the collective invasion model

Many model settings can be edited in the [PhysiCell_settings.xml](config/PhysiCell_settings.xml).  While you don't have to edit this file you can get different behavior by changing the settings in the file.  

You can follow these suggestions to familiarize yourself with the model.

1) Run this model using the default parameters. This will be produce the right hand side of the figure at the beginning of this README (Figure 5 of Reference [2]). 

2) Change anisotropy_increase_rate to 0.001 and fiber_realignment_rate to 1. This will decrease the collective behavior - eliminating collective invasion and leaving many followers in the center of the domain.

3) Change discrete_ECM_remodeling to 0. Note now the ECM is instantly remodeled - generating strong, clear signals for followers to read (view the anisotropy field). This will recover behavior similar to the default parameters and produce the left side of the figure included above.

4) Change default_cell_speed to 1.0 and then 0.25 - producing stigmergy and a non-changing morphology.

5) Now, you can enjoy yourself changing other parameters and creating new responses.

Note that code does not need recompiled in between parameter changes; the executable will parse the changes to the xml. 



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


## Acknowledgements

This work was funded in part by a joint (AMIGOS) JKTGF and BCRF grant. Thank you Margherita Botticelli for many productive conversations on the cell-ECM interaction code.

## References

[1] Metzcar, J, Duggan, BS, Fischer, B, Murphy, M, Heiland, R., Macklin, P. A simple framework for agent-based
modeling with extracellular matrix. bioRxiv 2022.11.21.514608; doi: [https://doi.org/10.1101/2022.11.21.514608](https://doi.org/10.1101/2022.11.21.514608)

[2] Ghaffarizadeh, A, Heiland, R, Friedman, SH, Mumenthaler, SM, and Macklin, P. PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: [10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)

[3] Cheung KJ, Gabrielson E, Werb Z, Ewald AJ. Collective invasion in breast cancer requires a conserved basal epithelial program. Cell 2013; 155(7):1639-51. [10.1016/j.cell.2013.11.029](10.1016/j.cell.2013.11.029)

[4] Nguyen-Ngoc KV\*, Cheung KJ*, Brenot A, Shamir ER, Gray RS, Hines WC, Yaswen P, Werb Z, Ewald AJ. The ECM microenvironment regulates collective migration and local dissemination in normal and malignant mammary epithelium. Proceedings of the National Academy of Science 2012; [10.1073/pnas.1212834109](10.1073/pnas.1212834109) *Co-First Authors. PMCID: PMC3465416

[5] Duggan, BS, Metzcar, J, and Macklin, P (2021). DAPT: A package enabling distributed automated parameter testing. Gigabyte 2021, 1–10. [10.46471/gigabyte.22](10.46471/gigabyte.22).


**Latest PhysiCell info:**  follow [@PhysiCell](https://twitter.com/PhysiCell) on Twitter (http://twitter.com/PhysiCell)
 


## Brief detailed description of ECM model and cell-ECM interactions


As mentioned above, the ECM model has three components: `anisotropy` (average local fiber-fiber alignment correlation), `fiber orientation` (average orientation of the fibers), and `density` (average relative volume fraction of ECM fibers). We place an array of ECM elements to spatially model the ECM in a tissue.

- `Density` (a scalar ranging from 0-1) volume fraction of fibers, which represents the average local  fiber density (range 0 - 1). Zero refers to completely fluid  filled and one to completely packed with  fibers without void space
- `Orientation` (numerically given as a unit vector) represents the overall (average) fiber orientation
- `Anisotropy` (0-1) average local  ber- ber alignment correlation (range 0 - 1). At zero, there is no correlation and at one, locally there is complete fiber-to-fiber correlation. 


_ECM impacts cell migration_:
- Fiber orientation provides directional cues.
- Anisotropy gives a strength of directional cue: high anisotropy increases an ECM element's influence on direction of cell migration.
- Density influences cell speed: too little ECM, cells have nothing to attach to; too much, cells cannot pass.

_Cell migration and movement impact microstructure_:
- Direction of cell migration reorients an ECM elements's orientation.
- Cell-ECM element contact increases ECM anisotropy proportional to cell speed.
- Cells remodel ECM density towards a target value.

This model is motivated by findings in the developmental, disease, and tumor biology literature as well as inspired by previous modeling e orts [1, 8, 9, 53, 61]. The cell-ECM interactions are specified at the cellular level, enabling a variety of cell-ECM interactions, in particular changes in cell motility and ECM remodeling capabilities [6, 9, 61]. Additionally, ECM variables can be used to impact other cellular behaviors such as proliferation and death. Finally, we note and will demonstrate that these features can be integrated with others such as sensitivity to chemical cues and cell-cell adhesion to obtain an even richer range of cell behaviors.




An `anisotropy` value closer to 0 means the fibers are less aligned, producing little signal, and values closer to 1 indicatd highly aligned fibers that produce a strong signal. `Density` influences cell speed as the cell gets stuck if the fibers are too thick and if the fibers are too sparse then the cells have nothing to "grab".  The `fiber orientation` influences the motility vector of the cells and the amount of influence is controlled by the `anisotropy`. ECM modifying cells can align the fibers in the direction of their movement (changing `fiber orientation`), increase the `anisotropy` in proportion to the cell's speed, and modify `density` up or down to a target density. 