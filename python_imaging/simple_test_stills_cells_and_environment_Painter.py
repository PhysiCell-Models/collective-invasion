import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

options_for_figure2d = {} # should be the same as 2b and 2c

options_for_figure2d = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : False,
                        "plot_ECM_anisotropy" : False,
                        "plot_ECM_orientation" : True,
                        "retrieve_ECM_data": True,
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : True,
                        "load_SVG_data": True,
                        'plot_cells_from_physicell_data': False,
                        'plot_cell_histogram': True,
                        'cell_alpha': 0.25,
                        "quiver_options" : {"scale_quiver": False,
                                          "mask_quiver": False},
                        "histogram_options" : {"num_bins": 40,
                                                  "vmax": 80,
                                                  "alpha_value": 1.0,
                                                  "color_bar": True,
                                                  "make_inset": True}
                       }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... I wonder if there
#### is a better/more robust way to do it (kwargs???, lots of "trapping"??) but this will be handled later ... and I can ask Randy etc
### What is up with scaling - hum ...


# if histogram_options is None:
#     num_bins = 40
#     vmax = 100
#     alpha_value = 1.0

mf = PhysiCellPlotter()

make_inset_flag = sys.argv[1] # just use as a string - its easier and more explicit

if make_inset_flag == 'make_inset':
    options_for_figure2d['histogram_options']['make_inset'] = True

elif make_inset_flag == 'dont_make_inset':
    options_for_figure2d['histogram_options']['make_inset'] = False

mf.generic_plotter(starting_index=600, number_of_samples=1, options=options_for_figure2d, file_name='still_contour_and_cells', )

# image_list_for_figure2d = []
# image_list_for_figure2d = [int(sys.argv[1]), int(sys.argv[2])]
# # image_list_for_figure2d = [150, 417]

# options_for_figure2d['plot_cell_histogram'] = False
# options_for_figure2d['plot_cells_from_SVG'] = True

# mf.generic_plotter(starting_index=699, number_of_samples=1, options=options_for_figure2d, file_name='still_cells')
# options_for_figure2d['load_full_physicell_data'] = False

# for number in image_list_for_figure2d: #### You'll have to update this in teh raw code!!!
#     mf.generic_plotter(starting_index=number, number_of_samples=number, options=options_for_figure2d, file_name='still_' + str(number))