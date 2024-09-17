# Python Visualizations

`image_processing_for_physicell.py` contains the class `PhysiCellPlotter`. This class enables reproducible visualization with many options - including setting which layers to include (cell, ECM components, diffusing fields), positional tracking of cell movement history, and smoothly generating stills and movies sequentially. It uses a modified version of [`pyMCDS.py`](https://github.com/PhysiCell-Tools/python-loader) - `pyMCDS_ECM.py`. The modifications include loading ECM data .mat outputs and formatting that data for use in visualizations. 

See the folder `scripts` for usage examples of the scripts. 

Note that due to continous development and iteration on the plots for the article accompanying this ECM extention, some now out of date figure numbers have crystallized into the the scripts included here. To the best of our ability, we include functional notes that include broadly which images were produced with which scripts. Future work includes smoothing these now historical names to broader/more generic names.  

## Random extra scripts

`run_movie_scripts_orthogonal.py` is in the folder `scripts` and is used to invoke the making of the movies for the stochastic replicates by invoking the appropriate video making script multiple times. This specific example is for the invasive cellular front orthgonal scenario, with similar scripts used to make the videos for the other stochastic replicates. 

Finally, there is a data vizualization in the paper accompanying this code. The code to produce that (Figure Supplementary Figure 5) is in `scripts`. It is `invasive_front_vizualization_of_variance.py` and uses output from `data_analysis_Painter.py`. 

## Simulation visualization scripts

Advanced scripts - using the PhysiCellPlotter class in `image_processing_for_physicell.py`. Place scripts in `output` or other directory at same level or alter path assignment. Examples used to produce stills and movies for [this preprint](https://www.biorxiv.org/content/10.1101/2022.11.21.514608) are given below. Using the various options available in the PhysiCellPlotter class, other visualizations can be created. 
* `full_history_movie_script.py`: Generates movie using all available outputs and plots the full cell positional histories as well as current cell positions and ECM fiber element orientations. 
full history still 
* `full_history_still_script.py`: Generates still, ploting the full cell positional histories up to the point of the still as well as current cell positions and ECM fiber element orientations
* `generic_movie_maker.py`: Generates a movie with default vizualization values. 
sys.path.append(r'../python_imaging')
*  `image_processing_script.py`: Generic example plotter. Outputs and shows a plot of the diffusing substrate and overlays cell positions. Loads frame 100. 
* `partial_history_2_level_contour_movie.py`: Generates of movie with each frame a recent portion of cell positional histories and anisotropy with 2 contour color levels (used for cases of instant ECM remodeling) as well as current cell positions and ECM fiber element orientations.
* `partial_history_2_level_contour_still.py`: Plots a recent portion of cell positional histories and anisotropy with 2 contour color levels (used for cases of instant ECM remodeling) as well as current cell positions and ECM fiber element orientations.
* `partial_history_multilevel_contour_movie.py`: Generates of movie with each frame a recent portion of cell positional histories and anisotropy with multiple contour color levels as well as current cell positions and ECM fiber element orientations.
* `partial_history_multilevel_contour_movie_density.py` : Deprecated - no longer used. 
* `partial_history_multilevel_contour_movie_fibrosis.py` : Generates movie for fibrosis simulation.
* `partial_history_multilevel_contour_movie_invasive_carcinoma.py` : Generates movie for basement membrane simulation. 
* `partial_history_multilevel_contour_still.py`: Plots a recent portion of cell positional histories and anisotropy with multiple contour color levels as well as current cell positions and ECM fiber element orientations.
* `partial_history_multilevel_contour_still_density.py` : Plots a recent portion of cell positional histories, a combination of anisotropy and ECM density with contour color levels specified at the command line as well as current cell positions and ECM fiber element orientations. Use in basement menbrane degradation stills. 
* `simple_test_march_movie.py` : Makes movie and series of stills that include anisotropy, element orientations, and cells. Used for cell march test. 
* `simple_test_movies_cells_and_environment_Painter.py` : Makes movie and series of still for the invasive cellular front simulations. 
* `simple_test_movies_cells_only.py` : Plots cells and cell positional histories, making every 10 cell red for visualization purposes. Prodcuces movie.
* `simple_test_stills_cells_and_environment_only.py` : Plots initial cell positions and element orientations and two additional stills (specified as output number at the command line) of cells and cell positional histories, making every 10 cell red for visualization purposes.
* `simple_test_stills_cells_and_environment_Painter.py` : Plots stills for the invasive cellular front experiment. Uses command line flag to plot inset of profile of cell positions. Used with `make_Painter_stills.sh` to produce plots for invase cellular front plot. 
* `simple_test_stills_march.py` : Plots the stills for the "cell march" figure - including element orientations, anisotropy, and still positions. 



Basic scripts. Original basis for integrated plotter in `image_processing_for_physicell.py`. Included as basic original examples. Place scripts in `output` or other directory at same level or alter path assignment. 
* `basic_cell_plot.py`: create a plot with Cells plotted. X
* `cell_march_movie_script.py` : Near copy of simple_test_march_movie. Deprecated. 
* `cell_march_stills_script.py` :  Near copy of simple_test_stills_march. Deprecated. 
* `cell_plus_environment_movie_maker.py`: Generates a movie of the ECM anisotorpy and orientations with cell overlay. Alter as needed for need. 
* `cell_plus_environment_plotter.py`: Generates a still of the ECM anisotorpy and orientations with cell overlay. Alter as needed for need. 
* `cell_track_plotter.py`: Plots still of cells and cell positional histories. Basis for cell positional history plotting function. 
* `cell_tracker_movie.py`: Generates a movie of cells and cell positional histories. 
* `finished-combined-plot.py`: Generates a layered plot combining cell positions and multiple aspcts of the microenvironment. 

## Future work 

Future work in python_imaging includes changes from pyMCDS based-parsing to [pcdl-based](https://github.com/PhysiCell-Tools/python-loader) data handling and reading. 

Remove variable names that reference specific, out of data figure numbers.

Remove hard coding of cell radii. 