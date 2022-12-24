import sys
import matplotlib.pyplot as plt

# To use the imaging module wihtout having to put it in every directory for analysis, put in the absolute path below. That directory
# will also need either 'pyMCDS.py' or 'pyMCDS_ECM.py'. Otherwise, place 'imaging_processing_for_phyiscell.py' and either MCDS file 
# in the current working directory. 

# sys.path.append(r'')

from image_processing_for_physicell import *

options_for_figure = {}

options_for_figure = {"output_plot" : True,
                       "show_plot" : True,
                       "produce_for_panel" : False,
                        "plot_ECM_anisotropy" : False,
                        "plot_ECM_orientation" : False,
                        "retrieve_ECM_data": False,
                        "retrieve_first_chemical_field_data" : True,
                        'plot_chemical_field' : True, 
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : True,
                        "contour_options" : {'lowest_contour': 0.0, ### I woud like this to be cleaner - but it does work!!!
                                           'upper_contour': 38,
                                           'number_of_levels': 38,
                                           'color_map_name': 'summer',
                                           'color_bar': True
                                           },
                        "quiver_options" : None
                       }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... I wonder if there
#### is a better/more robust way to do it (kwargs???, lots of "trapping"??) but this will be handled later ... and I can ask Randy etc
### What is up with scaling - hum ...

# oof - I got different results on the two runs when I did and didn't scale by anistoropy ... yikes! How do I manage that!!

mf = PhysiCellPlotter()

# options['']

# plot_cells_and_uE_for_movie (0, 1, 10, 1999)
#
#
# plot_cell_tracks_from_svg(0, 1, 10)
# general_image_plotter (filename: str, folder: str='.', output_folder='', cell_df: dict=None, cell_positions_from_SVG: dict=None, chemical_mesh: dict=None, ECM_mesh: dict=None, options=None):

image_list_for_figure = []

image_list_for_figure = [100]

for number in image_list_for_figure:
    mf.generic_plotter(starting_index=0, number_of_samples=number, options=options_for_figure)
# generic_plotter (start, intervnal, finish, save_filename, data_path, save_path, options)
#
#     All based on options/logic- function
#     load_cell_positiondata
#     load_uE_data_chemical
#     load_uE_data_ECM
#
#     process data into plots - functions
#     - cell tracks (might be loaded by just be plotted???)
#     - cell positions
#     - ECM layer
#     - chemical layer
#
#     complete plot presentaiont and save (maybe functions)
#     - title
#     - axes