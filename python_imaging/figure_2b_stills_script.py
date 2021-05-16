import sys
import matplotlib.pyplot as plt
sys.path.append(r'/Users/JohnMetzcar/Documents/GitHub/AMIGOS-invasion/python_imaging')

from image_processing_for_physicell import *

options_for_figure2b = {}

options_for_figure2b = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : True,
                        "plot_ECM_anisotropy" : False,
                        "plot_ECM_orientation" : True,
                        "retrieve_ECM_data": True,
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : True,
                        "load_SVG_data": True,
                        'plot_cells_from_physicell_data': False,
                        "quiver_options" : {"scale_quiver": False,
                                          "mask_quiver": False}
                       }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... I wonder if there
#### is a better/more robust way to do it (kwargs???, lots of "trapping"??) but this will be handled later ... and I can ask Randy etc
### What is up with scaling - hum ...

# oof - I got different results on the two runs when I did and didn't scale by anistoropy ... yikes! How do I manage that!!

mf = PhysiCellPlotter()
# m2 = PhysiCellPlotter()
# m3 = PhysiCellPlotter()

mf.generic_plotter(starting_index=0, number_of_samples=1, options=options_for_figure2b, file_name='horizontal_ECM_w_chemical_cue_0')

image_list_for_figure2b = []

image_list_for_figure2b = [50, 150]

# file_name = 'horizontal_ECM_w_chemical_cue_' + str(90)

options_for_figure2b['plot_ECM_orientation'] = False

for number in image_list_for_figure2b:
    mf.generic_plotter(starting_index=0, number_of_samples=number, options=options_for_figure2b, file_name='horizontal_ECM_w_chemical_cue_' + str(number))

# mf.generic_plotter(starting_index=90, number_of_samples=1, options=options_for_figure2a)
# m2.generic_plotter(starting_index=500, number_of_samples=1, options=options_for_figure2a)
# m3.generic_plotter(starting_index=1200, number_of_samples=1, options=options_for_figure2a)

# mf.generic_plotter (number_of_samples=10, options=options_for_figure2a)
# mf.create_separate_colorbar(contour_options=options_for_figure2a['contour_options'])

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