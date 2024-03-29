import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

options_for_figure2d = {} # should be the same as 2b and 2c

options_for_figure2d = {"output_plot" : True,
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

mf = PhysiCellPlotter()

mf.generic_plotter(starting_index=0, number_of_samples=1, options=options_for_figure2d, file_name='still_initial_step')

image_list_for_figure2d = []
image_list_for_figure2d = [int(sys.argv[1]), int(sys.argv[2])]
# image_list_for_figure2d = [150, 417]

options_for_figure2d['plot_ECM_orientation'] = False
options_for_figure2d['retrieve_ECM_data'] = False
options_for_figure2d['load_full_physicell_data'] = False

for number in image_list_for_figure2d:
    mf.generic_plotter(starting_index=0, number_of_samples=number, options=options_for_figure2d, file_name='still_' + str(number))