import sys
import matplotlib.pyplot as plt
sys.path.append(r'/Users/JohnMetzcar/Documents/GitHub/AMIGOS-invasion/python_imaging')

from image_processing_for_physicell import *

options_for_figure2a = {}

options_for_figure2a = {"output_plot" : True,
                       "show_plot" : True,
                       "produce_for_panel" : False,
                        "plot_ECM_anisotropy" : True,
                        "plot_ECM_orientation" : True,
                        "plot_cells_from_SVG" : True
                       }

mf = PhysiCellPlotter()

# options['']

# plot_cells_and_uE_for_movie (0, 1, 10, 1999)
#
#
# plot_cell_tracks_from_svg(0, 1, 10)
# general_image_plotter (filename: str, folder: str='.', output_folder='', cell_df: dict=None, cell_positions_from_SVG: dict=None, chemical_mesh: dict=None, ECM_mesh: dict=None, options=None):
mf.generic_plotter (number_of_samples=10, options=options_for_figure2a)

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