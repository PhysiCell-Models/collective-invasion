# Python Visualizations

`image_processing_for_physicell.py` contains the class `PhysiCellPlotter`. This class enables reproducible visualization with many options - including setting which layers to include (cell, ECM components, diffusing fields), positional tracking of cell movement history, and smoothly generating stills and movies sequentially. It uses a modified version of [`pyMCDS.py`](https://github.com/PhysiCell-Tools/python-loader) - `pyMCDS_ECM.py`. The modifications include loading ECM data .mat outputs and formatting that data for use in visualizations. 


## Scripts



Advanced scripts - using the PhysiCellPlotter class in `image_processing_for_physicell.py`. Place scripts in `output` or other directory at same level or alter path assignment. 
* 


Basic scripts. Original basis for integrated plotter in `image_processing_for_physicell.py`. Included as basic original examples. Place scripts in `output` or other directory at same level or alter path assignment. 
* `basic_cell_plot.py`: create a plot with Cells plotted.
* `cell_plus_environment_movie_maker.py`: Generates a movie of the ECM anisotorpy and orientations with cell overlay. Alter as needed for need. 
* `cell_plus_environment_plotter.py`: Generates a still of the ECM anisotorpy and orientations with cell overlay. Alter as needed for need. 
* `cell_track_plotter.py`: Plots still of cells and cell positional histories. 
* `cell_tracker_movie.py`: Generates a movie of cells and cell positional histories.
