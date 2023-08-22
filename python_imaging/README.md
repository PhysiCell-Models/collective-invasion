# Python Visualizations

`image_processing_for_physicell.py` contains the class `PhysiCellPlotter`. This class enables reproducible visualization with many options - including setting which layers to include (cell, ECM components, diffusing fields), positional tracking of cell movement history, and smoothly generating stills and movies sequentially. It uses a modified version of [`pyMCDS.py`](https://github.com/PhysiCell-Tools/python-loader) - `pyMCDS_ECM.py`. The modifications include loading ECM data .mat outputs and formatting that data for use in visualizations. 


## Scripts

Advanced scripts - using the PhysiCellPlotter class in `image_processing_for_physicell.py`. Place scripts in `output` or other directory at same level or alter path assignment. Examples used to produce stills and movies for [this preprint](https://www.biorxiv.org/content/10.1101/2022.11.21.514608) are given below. Using the various options available in the PhysiCellPlotter class, other visualizations can be created. 
* `full_history_movie_script.py`: Generates movie using all available outputs and plots the full cell positional histories as well as current cell positions and ECM fiber element orientations. 
* `generic_movie_maker.py`: Generates a movie with default vizualization values. 
sys.path.append(r'../python_imaging')
* `full_history_still_script.py`: Generates still, ploting the full cell positional histories up to the point of the still as well as current cell positions and ECM fiber element orientations
*  `image_processing_script.py`: Generic example plotter. Outputs and shows a plot of the diffusing substrate and overlays cell positions. Loads frame 100. 
* `partial_history_2_level_contour_movie.py`: Generates of movie with each frame a recent portion of cell positional histories and anisotropy with 2 contour color levels (used for cases of instant ECM remodeling) as well as current cell positions and ECM fiber element orientations.
* `partial_history_2_level_contour_still.py`: Plots a recent portion of cell positional histories and anisotropy with 2 contour color levels (used for cases of instant ECM remodeling) as well as current cell positions and ECM fiber element orientations.
* `partial_history_multilevel_contour_still.py`: Plots a recent portion of cell positional histories and anisotropy with multiple contour color levels as well as current cell positions and ECM fiber element orientations.
* `partial_history_multilevel_contour_movie.py`: Generates of movie with each frame a recent portion of cell positional histories and anisotropy with multiple contour color levels as well as current cell positions and ECM fiber element orientations.
* `simple_test_movies_cells_only.py` : Plots cells and cell positional histories, making every 10 cell red for visualization purposes. Prodcuces movie.
* `simple_test_stills_cells_only.py` : Plots cells and cell positional histories, making every 10 cell red for visualization purposes. 


Basic scripts. Original basis for integrated plotter in `image_processing_for_physicell.py`. Included as basic original examples. Place scripts in `output` or other directory at same level or alter path assignment. 
* `basic_cell_plot.py`: create a plot with Cells plotted.
* `cell_plus_environment_movie_maker.py`: Generates a movie of the ECM anisotorpy and orientations with cell overlay. Alter as needed for need. 
* `cell_plus_environment_plotter.py`: Generates a still of the ECM anisotorpy and orientations with cell overlay. Alter as needed for need. 
* `cell_track_plotter.py`: Plots still of cells and cell positional histories. 
* `cell_tracker_movie.py`: Generates a movie of cells and cell positional histories.
* `finished-combined-plot.py`: Generates a layered plot combining cell positions and multiple aspcts of the microenvironment. 
