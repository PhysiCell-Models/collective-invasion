import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

options_for_figure2d = {} # should be the same as 2b and 2c

options_for_figure2d = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : False,
                        "plot_ECM_anisotropy" : False,
                        "plot_ECM_orientation" : False,
                        "retrieve_ECM_data": False,
                        "load_full_physicell_data" : False,
                        "plot_cells_from_SVG" : True,
                        "load_SVG_data": True,
                        'plot_cells_from_physicell_data': False,
                        "produce_for_movie" : True,
                        "quiver_options" : {"scale_quiver": False,
                                          "mask_quiver": False}
                       }

movie_options_for_figure_2d = {}

movie_options_for_figure_2d = {'INCLUDE_ALL_SVGs': True,
                            'INCLUDE_FULL_HISTORY': True
                            }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... 

mf = PhysiCellPlotter()

mf.produce_movie(save_name='figure_2d_circular_ecm_with_chemotaxsis', movie_options=movie_options_for_figure_2d, image_options=options_for_figure2d)