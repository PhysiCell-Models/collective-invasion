import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

options_for_figure5c = {}

options_for_figure5c = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : True,
                        "plot_ECM_anisotropy" : True,
                        "plot_ECM_density" : False,
                        "plot_ECM_orientation" : True,
                        "retrieve_ECM_data": True,
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : True,
                        "load_SVG_data": True,
                        'plot_cells_from_physicell_data': False,
                        "contour_options" : {'lowest_contour': 1e-14, ### I woud like this to be cleaner - but it does work!!!
                                           'upper_contour': 1.0,
                                           'number_of_levels': 25,
                                           'color_map_name': 'Reds',
                                           'color_bar': False
                                           },
                       }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... 

mf = PhysiCellPlotter()

# in future, could iterate over input arguments to make this more general
image_list_for_figure5c = [int(sys.argv[1]), int(sys.argv[2])]

# 3rd argument is for plotting ECM density or anisotropy - but prior scripts that use ONLY 2 arguments will still work
if (len(sys.argv) > 3):
    if (sys.argv[3] == "plot_density"):
        options_for_figure5c["plot_ECM_density"] = True
        options_for_figure5c["plot_ECM_anisotropy"] = False
    # else, you get anisotropy
    
number_of_samples = 12
for number in image_list_for_figure5c:
    starting_index = number-number_of_samples + 1
    mf.generic_plotter(starting_index=starting_index, number_of_samples=number_of_samples, options=options_for_figure5c, file_name='multi_contour_still_' + str(number))
# trail_length must be at least 1!!!!!!!!
mf.create_separate_colorbar(contour_options = options_for_figure5c["contour_options"])

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