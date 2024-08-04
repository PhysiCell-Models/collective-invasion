import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

options_for_figure2a = {}

options_for_figure2a = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : True,
                        "plot_ECM_anisotropy" : True,
                        "plot_ECM_orientation" : True,
                        "retrieve_ECM_data": True,
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : False,
                        "load_SVG_data": False,
                        'plot_cells_from_physicell_data': True,
                        "contour_options" : {'lowest_contour': 0.90, ### I woud like this to be cleaner - but it does work!!!
                                           'upper_contour': 0.92,
                                           'number_of_levels': 25,
                                           'color_map_name': 'Reds',
                                           'color_bar': True
                                           },
                        "quiver_options" : {"scale_quiver": False,
                                          "mask_quiver": False}
                       }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... I wonder if there
#### is a better/more robust way to do it (kwargs???, lots of "trapping"??) but this will be handled later ... 

mf = PhysiCellPlotter()

image_list_for_figure2a = []

image_list_for_figure2a = [90, 500, 1200]

file_name = 'march_' + str(90)
print('To make \'March\' images - must hard code both colors to blue!!!!')
print('SEE plot_cells_from_physicell_data function!!!')
for number in image_list_for_figure2a:
    mf.generic_plotter(starting_index=number, number_of_samples=1, options=options_for_figure2a, file_name='march_' + str(number))

