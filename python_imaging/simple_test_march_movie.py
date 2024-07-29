import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

options_for_figure2a = {}

options_for_figure2a = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : False,
                        "plot_ECM_anisotropy" : True,
                        "plot_ECM_orientation" : True,
                        "retrieve_ECM_data": True,
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : False,
                        "load_SVG_data": False,
                        'plot_cells_from_physicell_data': True,
                        "produce_for_movie" : True,
                        "contour_options" : {'lowest_contour': 0.90, ### I woud like this to be cleaner - but it does work!!!
                                           'upper_contour': 0.92,
                                           'number_of_levels': 25,
                                           'color_map_name': 'Reds',
                                           'color_bar': True
                                           },
                        "quiver_options" : {"scale_quiver": False,
                                          "mask_quiver": False}
                       }

movie_options_for_figure_2a = {}

movie_options_for_figure_2a = {'INCLUDE_ALL_SVGs': True,
                            'INCLUDE_FULL_HISTORY': False
                            }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... I wonder if there
#### is a better/more robust way to do it (kwargs???, lots of "trapping"??) but this will be handled later. 

mf = PhysiCellPlotter()

mf.produce_movie(save_name='simple_test_cells', movie_options=movie_options_for_figure_2a, image_options=options_for_figure2a)