import sys
import argparse
# defined command line options
# this also generates --help and error handling

import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--snapshot_IDs",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=int,
  default=[15],  # default if nothing is provided
)

# CLI=argparse.ArgumentParser()
# CLI.add_argument(
#   "plot_orientations",  # name on the CLI - drop the `--` for positional/required parameters
# #   nargs="*",  # 0 or more values expected => creates a list
#   type=bool,
#   default=True  # default if nothing is provided
# )


CLI.add_argument(
  "--contour_range",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=float,
  default=[0.0, 1.0],  # default if nothing is provided
)

# parse the command line
args = CLI.parse_args()
# access CLI options
print(args.snapshot_IDs)
print(args.contour_range)

options_for_figure5c = {}

options_for_figure5c = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : True,
                        "plot_ECM_anisotropy" : True,
                        "plot_ECM_density" : True,
                        "plot_ECM_orientation" : True,
                        "retrieve_ECM_data": True,
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : True,
                        "load_SVG_data": True,
                        'plot_cells_from_physicell_data': False,
                        "contour_options" : {'lowest_contour': 0.5, ### I woud like this to be cleaner - but it does work!!!
                                           'upper_contour': 1.0,
                                           'number_of_levels': 25,
                                           'color_map_name': 'YlOrRd',
                                           'color_bar': False,
                                           'alpha': 0.5
                                           },
                       }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting ... 

mf = PhysiCellPlotter()

# in future, could iterate over input arguments to make this more general
# image_list_for_figure5c = [args.snapshot_IDs]

# 3rd argument is for plotting ECM density or anisotropy - but prior scripts that use ONLY 2 arguments will still work
# if (len(sys.argv) > 3):
#     if (sys.argv[3] == "plot_density"):
#         options_for_figure5c["plot_ECM_density"] = True
#         options_for_figure5c["plot_ECM_anisotropy"] = False
    # else, you get anisotropy

# if (len(sys.argv) > 3):
options_for_figure5c["contour_options"]["lowest_contour"] = args.contour_range[0]
options_for_figure5c["contour_options"]["highest_contour"] = args.contour_range[1]

# options_for_figure5c["plot_ECM_orientation"] = args.plot_orientations
    
trail_length = 8
for snapshot_ID in args.snapshot_IDs:
    if (snapshot_ID - trail_length + 1 < 0):
        starting_index = 0
        modified_trail_length = snapshot_ID + 1 
        print("WARNING: trail length is too long for this snapshot ID")
        print(snapshot_ID)
        print(trail_length)
    else:
        starting_index = snapshot_ID-trail_length + 1
        modified_trail_length = trail_length
    mf.generic_plotter(starting_index=starting_index, number_of_samples=modified_trail_length, options=options_for_figure5c, file_name='multi_contour_still_' + str(snapshot_ID))
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