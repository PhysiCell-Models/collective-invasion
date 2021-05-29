import math, os, sys, re
import xml.etree.ElementTree as ET
import numpy as np

import matplotlib.pyplot as plt

######## If using on remote system, uncomment this line below to load correct matplotlib backend ################
# matplotlib.use('Agg')

import matplotlib.colors as mplc
import matplotlib.colorbar as colorbar
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle

import distutils.util

# from pyMCDS_ECM import *
try:
    from pyMCDS_ECM import *
except ImportError:
    from pyMCDS import *

# Features for PhysiImage module
# Plots cells, tracks (as vectors???), and at least one microenvironment feature (ECM or otherwise)
# Allows for fine grain control of rate of plotting of tracks - start, end and interval
# Allows for fine grain control of outputs - quality, for insets, for videos
# has scale bar (ideally)
# preserves correctly scaled cell diamteers - DONE! working with SVG loader if cells are constant size. Must use other one otherwise.
# preserves cell colors (only in SVGs!!!!!!!!) and also allows for that to be overridden if needed
# Gets mat size and time/slide number from images
# Allows you to specify a title and add time/slide number to it
# plots color bar that isn't stupidly large (can I get something that returns a figure and then lets me change it in a script????) Ormaybe I can just write a bunch of different ones or use flags. Something to make it easier than it currently is - which it is currently assine.
# Be able to specify an output directory (might want to check that it is exsists (or not - does Python give an error?))
# Add in module catch that says - ECM functionality will fail - load pyMCDS_ECM to use with ECM, otherwise your are fine

class PhysiCellPlotter():

    # # https://realpython.com/documenting-python-code/
    # https://stackoverflow.com/questions/37019744/is-there-a-consensus-what-should-be-documented-in-the-classes-and-init-docst
    
    def __init__(self, parent = None):

        """
        Initializes a plot using matplotlib pyplot subplot routine, returning a figure and axis handle to self. Provides a default figure size, title, and all 
        default options (self.default_options) required to make a plot that displays ONLY cell positions (and tracks if the generic_plotter is called with the appropriate variables).
        self.default_options is used by generic_plotter to fill in any option values not set when pass the options to generic_plotter

        """
        self.figsize_width_svg = 7.0
        self.figsize_height_svg = 7.0
        self.title = "title"
        self.fig, self.ax = plt.subplots(figsize=(self.figsize_width_svg, self.figsize_height_svg))
        self.default_options = {"output_plot": True,
                       "show_plot": True,
                       "produce_for_panel": False,
                       "load_SVG_data" : True, # cell color and positions
                       "load_full_physicell_data" : False, # The runs py_MCDS_ECM (ECM could be split out later if pyMCDS changes??)
                       "retrieve_first_chemical_field_data" : False, # Gets first chemical field from pyMCDS object. Eventually will probably want multiple sets of options - like "load this field" etc - maybe need an options class??
                       "retrieve_ECM_data": False, # Gets ECM data from pyMCDS object
                       "plot_ECM_anisotropy" : False, # Calls contour plotter with anisotropy as input
                        'plot_chemical_field': False, # Calls contour plotter with chemical field as input
                       "plot_ECM_orientation" : False, # calls quiver plotter with orientation as input
                       "plot_cells_from_SVG" : True, # plots cell positions and colors using data from SVGs
                       "plot_cells_from_physicell_data": False, # plots cell positions from pyMCDS --> will need more options if I want to specify colors ... currently set up to read color from SVG data
                       ####### Cell tracks are always plotted when calling plot_cells_from_svg - to not plot tracks - make the number of samples = 1 ...
                        "produce_for_movie" : False,
                        "contour_options": None,
                        "quiver_options": None}

    def generic_plotter(self, starting_index: int = 0, sample_step_interval: int = 1, number_of_samples: int = 120,
                        file_name: str = None, input_path: str= '.', output_path: str= '', naming_index: int=0, options=None):

        """
        Produces multlilayer image: allows for one cell layer, a contour layer (with colorbar), vector field, 
        and cell positional history, plotted as arrows (quiver plot) with final cell positions plotted as a cirle.
        Options passed through a dictionary (see class consctructor for example). 

        sample_step_interval * number_of_samples - starting_index yields the trail length in time steps. number_of_samples provides
        the number of intervals plotted per image. 

        Example: starting_index of 0, sample intervale of 1, and number of samples of 120 will produce a cell track 120 steps long, sampled at whatever rate the SVGs were produced, starting at 
        snapshot 0 going until snapshot 119. 

        Parameters
        ----------
        starting_index :
            Integer index of the PhysiCell SVG output to begin trackign at. Default is 0.
        sample_step_interval :
            Interval (number of time steps (SVGs)) to sample at. A value of 2 would add a tracking point for every other SVG. Default is 1.
        number_of_samples :
            Total Number of SVGs to process. Length of cell positional history. Number_of_samples * sample_step_interval provides the index of the final SVG to process. Default is 120.
        file_name : 
            Use to specify a non-default image output name. "produce_for_movie" option=True overrides both the default and given (if given) file name to allow for 
            required image names to make movie. Default is None, producing the default naming scheme. Example: for the default arguements: 0_1_120 (starting index, sample interval, number of samples). 
        input_path : 
            Sets input directory for .mat, xml, and SVG files. All data assumed to be in the same directory. Default values is the current/working directory.
            NOT CURRENTLY IMMPLEMENTED FOR SVGs!!!!!!!!!! In future versions, plan to use os.chdir, but want to set up logic to help with this.
        output_path :
            Sets image output location. Default is current/working directory. 
        naming_index : 
            Special use variable to specify expected and ordered file names required to make movie from multiple output images. Default is 0.
        options : 
            Diectinoary containing all options required to specify image to be produced. Default is None. Since the dictionary is requied, the default trigeers copying of the default_options, 
            specified in the PhysiCellPlotter default constructor. Basically, the defaults make an image with cells and cell histories only plotted. 
        
        Returns
        -------
        Null :
            Produces a png image using specified PhysiCell inputs etc as specified in the options dictionary. 

        """

        self.fig, self.ax = plt.subplots(figsize=(self.figsize_width_svg, self.figsize_height_svg))

        if options is None:
            options = {"output_plot": True,
                       "show_plot": True,
                       "produce_for_panel": False,
                       "load_SVG_data" : True, # cell color and positions
                       "load_full_physicell_data" : False, # The runs py_MCDS_ECM (ECM could be split out later if pyMCDS changes??)
                       "retrieve_first_chemical_field_data" : False, # Gets first chemical field from pyMCDS object. Eventually will probably want multiple sets of options - like "load this field" etc - maybe need an options class??
                       "retrieve_ECM_data": False, # Gets ECM data from pyMCDS object
                       "plot_ECM_anisotropy" : False, # Calls contour plotter with anisotropy as input
                       'plot_chemical_field' : False,
                       "plot_ECM_orientation" : False, # calls quiver plotter with orientation as input
                       "plot_cells_from_SVG" : True, # plots cell positions and colors using data from SVGs
                       "plot_cells_from_physicell_data": False, # plots cell positions from pyMCDS --> will need more options if I want to specify colors ... currently set up to read color from SVG data
                       ####### Cell tracks are always plotted when calling plot_cells_from_svg - to not plot tracks - make the number of samples = 1 ...
                       "produce_for_movie" : False,
                       "contour_options": None,
                       "quiver_options": None}
        else:
            for key in self.default_options.keys():
                if key in options.keys():
                    pass
                else: 
                    options[key] = self.default_options[key]
                    print(options[key]) ##### Add in something saying that defaults were used for this key value???. Then is there someway to get it to only do that once per call???
                    print(key)

        # print("Current Working Directory " , os.getcwd())
        # os.chdir("/home/varun/temp")

        if options["load_SVG_data"] is True:
            cell_positions, cell_attributes, title_str, plot_x_extend, plot_y_extend = self.load_cell_positions_from_SVG(
            starting_index, sample_step_interval, number_of_samples)
            print('Stil need to get input_path for SVGs working!!!')

        if options["load_SVG_data"] is False:
            endpoint = starting_index + sample_step_interval * number_of_samples - 1
            final_snapshot_name = 'output' + f'{endpoint:08}'
            print(final_snapshot_name)
            title_str = 'some one should add extracting the file name from teh .mat files or similar to the code!!!'
            plot_x_extend = 1000
            plot_y_extend = 1000
            print("WARNING!!!!!!!!!!! Plot extent is not dynamic!!!!!!!!!!!!!! Load from SVG to get dynamic changes OR change pyMCDS to get bounding box then change Load Physicell Data method!!!!!")

        else:
            endpoint = starting_index + sample_step_interval * number_of_samples - 1
            final_snapshot_name = 'output' + f'{endpoint:08}'
            print(final_snapshot_name)
            title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str

        if file_name is None:
            file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)

        if options["produce_for_movie"] is True: 
            file_name = snapshot = 'output' + f'{naming_index:08}'
            print('Output file name forced to indexable name to produce movie')

        if options['load_full_physicell_data'] is True:
            self.load_full_physicell_data(final_snapshot_name, folder=input_path)
            print('test input path option!!!! (for loading physicell data...)')

        if options['retrieve_first_chemical_field_data'] is True:
            xx, yy, plane_oxy = self.load_chemical_field('oxygen')
            print('this call needs updated to use an option for putting in the chemical field name then defaulting ot oxygen perhaps for generic???')
        if options['retrieve_ECM_data'] is True:
            xx_ecm, yy_ecm, ECM_anisotropy, ECM_density, ECM_x_orientation, ECM_y_orientation = self.retreive_ECM_data()
        # 1e-14, 1.0

        # contour_spacing = np.linspace(contour_options['lowest_contour'], contour_options['upper_contour'],
        #                               contour_options['number_of_levels'])
        #
        # cs = self.ax.contourf(x_mesh, y_mesh, data_to_contour, cmap=contour_options['color_map_name'],
        #                       levels=contour_spacing)

        # if contour_options['color_bar'] is True:
        if options['plot_chemical_field'] is True:
            self.create_contour_plot(x_mesh=xx, y_mesh=yy, data_to_contour=plane_oxy,
                                     contour_options=options["contour_options"], options=options)
        if options['plot_ECM_anisotropy'] is True:
            self.create_contour_plot(x_mesh=xx_ecm, y_mesh=yy_ecm, data_to_contour=ECM_anisotropy, contour_options=options["contour_options"], options=options)

        if options['plot_ECM_orientation'] is True:
            self.create_quiver_plot(scaling_values=ECM_anisotropy, x_mesh=xx_ecm, y_mesh=yy_ecm, x_orientation=ECM_x_orientation, y_orientation=ECM_y_orientation, quiver_options=options['quiver_options'])
            # Would be greato to pass kwargs here to teh plotting function, but can do that later ... I think maybe I can do some default behavior here??
            # And have a scaling inconsistency - but can deal with that later ...
            # https://stackoverflow.com/questions/49887526/rescaling-quiver-arrows-in-physical-units-consistent-to-the-aspect-ratio-of-the/49891134
        if options['plot_cells_from_physicell_data'] is True:
            self.plot_cells_from_physicell_data()

        if options['plot_cells_from_SVG'] is True:
            self.create_cell_layer_from_SVG(cell_positions, cell_attributes)

        self.plot_figure(title_str, plot_x_extend, plot_y_extend, file_name, output_path, options)

    def plot_cells_from_physicell_data(self):
        cell_df = self.mcds.get_cell_df()
        cell_df['radius'] = (cell_df['total_volume'].values * 3 / (4 * np.pi)) ** (1 / 3)
        types = cell_df['cell_type'].unique()
        colors = ['yellow', 'blue']
        print('WARNING!!!!!! WARNING!!!!!!!!!! These colors are hard coded AND WONT WORK ON NON-ECM SIMS!!!!!')
        # Add cells layer
        for i, ct in enumerate(types):
            plot_df = cell_df[cell_df['cell_type'] == ct]
            for j in plot_df.index:
                circ = Circle((plot_df.loc[j, 'position_x'], plot_df.loc[j, 'position_y']),
                              radius=plot_df.loc[j, 'radius'], color=colors[i], alpha=0.7)
                # for a blue circle with a black edge
                # circ = Circle((plot_df.loc[j, 'position_x'], plot_df.loc[j, 'position_y']),
                #                radius=plot_df.loc[j, 'radius'], alpha=0.7, edgecolor='black')
                self.ax.add_artist(circ)

    def load_full_physicell_data (self, snapshot: str='output000000000', folder: str='.'):
        # load cell and microenvironment data
        self.mcds = pyMCDS(snapshot + '.xml', folder)

        # loads and reads ECM data
        self.mcds.load_ecm(snapshot + '_ECM.mat', folder)

    def create_contour_plot(self, x_mesh: dict, y_mesh: dict, data_to_contour: dict, contour_options=None, options: dict=None):
        ### best options are probably to just allow defaults, search for max and min for limits, or maybe insist on limits ...
        ### another obvious option - and this coudl be a global to reset ... you could even change it with function calls
        ### countour color maps ...

        if contour_options is None:
            cs = self.ax.contourf(x_mesh, y_mesh, data_to_contour, cmap="Reds")
            self.fig.colorbar(cs, ax=self.ax)
            # self.fig.show()
        else:

            # Make levels for contours
            contour_spacing = np.linspace(contour_options['lowest_contour'], contour_options['upper_contour'], contour_options['number_of_levels'])

            cs = self.ax.contourf(x_mesh, y_mesh, data_to_contour, cmap=contour_options['color_map_name'], levels=contour_spacing)

            if contour_options['color_bar'] is True:
                divider = make_axes_locatable(self.ax)
                cax = divider.append_axes("right", size="5%", pad=0.10)
                # other fancy things you can do with colorbars - https://stackoverflow.com/questions/16595138/standalone-colorbar-matplotlib
                if options is None:
                    cb = self.fig.colorbar(cs, cax=cax, format='%.3f')
                elif options['produce_for_panel'] is False:
                    cb = self.fig.colorbar(cs, cax=cax, format='%.3f')
                else:
                    tick_spacing = np.linspace(contour_options['lowest_contour'], contour_options['upper_contour'], 5)
                    cb = self.fig.colorbar(cs, cax=cax, format='%.2f', ticks=tick_spacing)
                    cb.ax.tick_params(labelsize=20)

    def create_separate_colorbar(self, file_name='just_colorbar', contour_options: dict=None):
        print('Working - gives continous colorbar instead of discrete - could fix possibly but not sure how to match N')

        if contour_options is not None:
            contour_spacing = np.linspace(contour_options['lowest_contour'], contour_options['upper_contour'],
                                          contour_options['number_of_levels'])
            fig, ax = plt.subplots(figsize=(1, 8))
            cmap_str = 'mpl.cm.' + contour_options['color_map_name']

            cmap = eval(cmap_str)
            norm = mpl.colors.Normalize(vmin=contour_options['lowest_contour'], vmax=contour_options['upper_contour'])
            cb = colorbar.ColorbarBase(ax, orientation='vertical',
                                           cmap=cmap, norm=norm)

            plt.savefig(output_folder + file_name, bbox_inches='tight', dpi=256)
            plt.show()
        else:
            print("you need to put in something for the color bar options. Supply \"contour_options\" to me!!!!")

    def create_quiver_plot(self, scaling_values: dict, x_mesh: dict, y_mesh: dict, x_orientation: dict, y_orientation: dict, quiver_options: dict=None):
        
        if quiver_options is None:
            mask = scaling_values > 0.0001
            ECM_x = np.multiply(x_orientation, scaling_values)
            ECM_y = np.multiply(y_orientation, scaling_values)
            self.ax.quiver(x_mesh[mask], y_mesh[mask], ECM_x[mask], ECM_y[mask],
                           pivot='middle', angles='xy', scale_units='inches', scale=4.75, headwidth=0, alpha = 0.3)
        else:
            if quiver_options["scale_quiver"] is True:
                ECM_x = np.multiply(x_orientation, scaling_values)
                ECM_y = np.multiply(y_orientation, scaling_values)
            else:
                ECM_x = x_orientation
                ECM_y = y_orientation

            # q = ax.quiver(xx_ecm[mask], yy_ecm[mask], scaled_ECM_x[mask], scaled_ECM_y[mask], pivot='middle', angles='xy', scale_units='inches', scale=2.0, headwidth=0,
            #               width=0.0015)  ## What is the deal with the line segment lengths shifting as the plots progress when I don't ue teh scaling??

            # mask out zero vectors
            mask = scaling_values > 0.0001
            if quiver_options["mask_quiver"] is True:
                self.ax.quiver(x_mesh[mask], y_mesh[mask], ECM_x[mask], ECM_y[mask],
                               pivot='middle', angles='xy', scale_units='inches', scale=4.75, headwidth=0, alpha = 0.3)
            else:
                self.ax.quiver(x_mesh, y_mesh, ECM_x, ECM_y,
                pivot='middle', angles='xy', scale_units='inches', scale=4.75, headwidth=0, alpha = 0.3)

    def load_chemical_field(self, field_name: str=None):

        #### Diffusion microenvironment
        xx, yy = self.mcds.get_2D_mesh()  # Mesh

        if field_name is not None:
            scalar_field_at_z_equals_zero = self.mcds.get_concentrations(field_name, 0.0)  # Oxyen (used for contour plot)
        else:
            print('Must supply field name as a string to use \'load_chemical_field\' function.')


        return xx, yy, scalar_field_at_z_equals_zero

    def retreive_ECM_data(self):

        #### ECM microenvironment
        xx_ecm, yy_ecm = self.mcds.get_2D_ECM_mesh()  # Mesh
        anisotropy_at_z_equals_zero = self.mcds.get_ECM_field('anisotropy', 0.0)  # Anistropy (used for scaling and contour plot)
        density_at_z_equals_zero = self.mcds.get_ECM_field('density', 0.0)
        x_orientation_at_z_equals_zero = self.mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0]
        y_orientation_at_z_equals_zero = self.mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0]

        return xx_ecm, yy_ecm, anisotropy_at_z_equals_zero, density_at_z_equals_zero, x_orientation_at_z_equals_zero, y_orientation_at_z_equals_zero

    def plot_figure(self, title_str: str, plot_x_extend: float, plot_y_extend: float, file_name: str, output_directory: str='', options: dict=None): ###### This should probably have to have options??????? Why though???
        if options is None:
            options= {"output_plot": True,
                      "show_plot": True,
                      "produce_for_panel": False
                      }
        output_plot = options['output_plot']
        show_plot = options['show_plot']
        produce_for_panel = options['produce_for_panel']
        output_folder = ''
        # print(output_folder.type)
        # print(file_name.type)
        # fig.figure(figsize=(7, 7))

        self.ax.set_aspect("equal")
        # endpoint = starting_index + sample_step_interval*number_of_samples
        #### Build plot frame, titles, and save data
        self.ax.set_ylim(-plot_y_extend/2, plot_y_extend/2)
        self.ax.set_xlim(-plot_x_extend/2, plot_x_extend/2)
    
        if produce_for_panel == False:
            # title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str
            # %"Starting at frame {}, sample interval of {} for {} total samples".format(number_of_samples, sample_step_interval, number_of_samples)
            self.ax.set_title(title_str)
        else:
            self.ax.xaxis.set_tick_params(labelsize=20)
            self.ax.yaxis.set_tick_params(labelsize=20)
            self.ax.set_xlabel('microns', fontsize=20)
            self.ax.set_ylabel('microns', fontsize=20)
            self.ax.set_xticks([ -plot_x_extend/2, -plot_x_extend/4, 0, plot_x_extend/4 ,plot_x_extend/2])
            self.ax.set_yticks([ -plot_y_extend/2, -plot_y_extend/4, 0, plot_y_extend/4 ,plot_y_extend/2])
            self.fig.tight_layout()
        # could change to the custom in the movie output or some other more better output if desired.
        output_folder = output_directory
        # if file_name is None:
        #     file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
    
        # Produce plot following the available options.
    
        if output_plot is True:
            plt.savefig(output_folder + file_name + '.png', dpi=256)
        if show_plot is True:
            plt.show()
        # self.fig.clf()
    
    def load_cell_positions_from_SVG(self, starting_index: int, sample_step_interval: int, number_of_samples: int):
            """
            Produces savable image of cell positional history, plotted as arrows (quiver plot) with final cell positions plotted as a cirle.
            Slight modification of the function in cell_track_plotter. The modification allows for tracking the index of a series
            of inputs such that outputs of this function can be appropriate indexed and compiled into a movie.
    
            sample_step_interval * number_of_samples - starting_index yields the trail length in time steps. number_of_samples provides
            the number of intervals plotted per image.
    
            Parameters
            ----------
            starting_index :
                Integer index of the PhysiCell SVG output to begin trackign at
            sample_step_interval :
                Interval (number of time steps (SVGs)) to sample at. A value of 2 would add a tracking point for every other SVG
            number_of_samples :
                Number of SVGs to process (total)/Length of cell positional history. Number of samples * sample size step interval provides the index of the final SVG to process
            output_plot :
                Save plot flag (required to produce a movie from resulting images)
            show_plot :
                Show plot flag (for processing many images in a loop, this should likely be set to false. Images have to be closed manually)
            produce_for_panel :
                    Flag - calls tight_layout, increases axes font sizes, and plots without title. For using in panels of images where there will be captions.
            Returns
            -------
            Null :
                Produces a png image from the input PhysiCell SVGs.
    
            """
    
            # if options is None:
            #     options = {"output_plot": True,
            #                "show_plot": True,
            #                "produce_for_panel": False
            #                }
            # output_plot = options['output_plot']
            # show_plot = options['show_plot']
            # produce_for_panel = options['produce_for_panel']
    
            d = {}  # dictionary to hold all (x,y) positions of cells
            d_attributes = {}  # dictionary to hold other attributes, like color (a data frame might be nice here in the long run ... ) \
            # currently only being read once only as cell dictionary is populated - so only use for static values!
    
            """ 
            --- for example ---
            In [141]: d['cell1599'][0:3]
            Out[141]: 
            array([[ 4900.  ,  4900.  ],
                   [ 4934.17,  4487.91],
                   [ 4960.75,  4148.02]])
            """
    
            ####################################################################################################################
            ####################################            Generate list of file indices to load       ########################
            ####################################################################################################################
    
            endpoint = starting_index + sample_step_interval * number_of_samples
            file_indices = np.linspace(starting_index, endpoint, num=number_of_samples, endpoint=False)
            print(file_indices)
    
            maxCount = starting_index
    
            ####### Uncomment for statement below to generate a random list of file names versus the prespecifed list. ########
            ####### Leaving for historical record. If used, the inputs would need to be a single integer,       ########
            ####### versus the three integers required to generate the prespecified list. Also, remove the other for statement. ########
            # count = 0
            #
            # for fname in glob.glob('snapshot*.svg'):
            #     print(fname)
            # # for fname in['snapshot00000000.svg', 'snapshot00000001.svg']:
            # # for fname in['snapshot00000000.svg']:
            # #   print(fname)
            #   count += 1
            #   if count > maxCount:
            #     break
    
            ####################################################################################################################
            ####################################        Main loading and processing loop                ########################
            ####################################################################################################################
    
            for file_index in file_indices:
                print(os.getcwd())
                fname = "%0.8d" % file_index
                fname = 'snapshot' + fname + '.svg'  # https://realpython.com/python-f-strings/
                print(fname)
    
                ##### Parse XML tree into a dictionary called 'tree" and get root
                # print('\n---- ' + fname + ':')
                tree = ET.parse(fname)
    
                # print('--- root.tag, root.attrib ---')
                root = tree.getroot()
                #  print('--- root.tag ---')
                #  print(root.tag)
                #  print('--- root.attrib ---')
                #  print(root.attrib)
                #  print('--- child.tag, child.attrib ---')
    
                numChildren = 0
    
                ### Find branches coming from root - tissue parents
                for child in root:
                    # print(child.tag, child.attrib)
                    #    print('attrib=',child.attrib)
                    #  if (child.attrib['id'] == 'tissue'):
    
                    if child.text and "Current time" in child.text:
                        svals = child.text.split()
                        title_str = "Current time: " + svals[2] + "d, " + svals[4] + "h, " + svals[
                            7] + "m"
    
                    if 'width' in child.attrib.keys():
                        #### Assumes 100% of difference in SVG width and height is due to top margin of the SVG!!!!!!
                        # print('Reading SVG - Assumes 100% of difference in SVG width and height is due to top margin of the SVG!!!!!!')
                        plot_x_extend = float(child.attrib['width'])
                        top_margin_size = abs(float(child.attrib['height']) - float(child.attrib['width']))

                        #### Remove the padding placed into the SVG to determine the true y extend
                        plot_y_extend = float(child.attrib['height']) - top_margin_size
    
                        #### Find the coordinate transform amounts
                        y_coordinate_transform = plot_y_extend / 2
                        x_coordinate_transform = plot_x_extend / 2
    
                    ##### Find the tissue tag and make it child
                    if 'id' in child.attrib.keys():
                        #      print('-------- found tissue!!')
                        tissue_parent = child
                        break
    
                #  print('------ search tissue')
    
                ### find the branch with the cells "id=cells" among all the branches in the XML root
                for child in tissue_parent:
                    #    print('attrib=',child.attrib)
                    if (child.attrib['id'] == 'cells'):
                        #      print('-------- found cells!!')
                        cells_parent = child
                        break
                    numChildren += 1
    
                ### Search within the cells branch for all indiviual cells. Get their locations
                num_cells = 0
                #  print('------ search cells')
                for child in cells_parent:
                    #    print(child.tag, child.attrib)
                    #    print('attrib=',child.attrib)
    
                    # Find the locations of the cells within the cell tags
                    for circle in child:
                        #      print('  --- cx,cy=',circle.attrib['cx'],circle.attrib['cy'])
                        xval = float(circle.attrib['cx'])
    
                        # should we test for bogus x,y locations??
                        if (math.fabs(xval) > 10000.):
                            print("xval=", xval)
                            break
                        yval = float(circle.attrib['cy'])  # - y_coordinate_transform
                        if (math.fabs(yval) > 10000.):
                            print("yval=", yval)
                            break
    
                        # Pull out the cell's location. If ID not already in stack to track, put in new cell in dictionary while applying coordinate transform.
                        if (child.attrib['id'] in d.keys()):
                            d[child.attrib['id']] = np.vstack((d[child.attrib['id']],
                                                               [float(circle.attrib['cx']) - x_coordinate_transform,
                                                                float(circle.attrib['cy']) - y_coordinate_transform]))
                        #### Comment out this else to produce single cell tracks
                        else:
                            d[child.attrib['id']] = np.array([float(circle.attrib['cx']) - x_coordinate_transform,
                                                              float(circle.attrib['cy']) - y_coordinate_transform])
                            d_attributes[child.attrib['id']] = circle.attrib['fill']
    
                        ###### Uncomment this elif and else to produce single cell tracks
                        # elif (child.attrib['id'] == 'cell24'):
                        #     d[child.attrib['id']] = np.array( [ float(circle.attrib['cx'])-x_coordinate_transform, float(circle.attrib['cy'])-y_coordinate_transform])
                        #     d_attributes[child.attrib['id']] = circle.attrib['fill']
                        # else:
                        #     break
    
                        ##### This 'break' statement is required to skip the nucleus circle. There are two circle attributes. \
                        ##### If both nuclear and cell boundary attributes are needed, this break NEEDS REMOVED!!!!
                        break
    
                        ### Code to translate string based coloring to rgb coloring. Use as needed.
                        # s = circle.attrib['fill']
                        # print("s=",s)
                        # print("type(s)=",type(s))
                        # if (s[0:3] == "rgb"):  # if an rgb string, e.g. "rgb(175,175,80)"
                        #     #  circle.attrib={'cx': '1085.59','cy': '1225.24','fill': 'rgb(159,159,96)','r': '6.67717','stroke': 'rgb(159,159,96)','stroke-width': '0.5'}
                        #     rgb = list(map(int, s[4:-1].split(",")))
                        #     rgb[:] = [x / 255. for x in rgb]
                        # else:  # otherwise, must be a color name
                        #     rgb_tuple = mplc.to_rgb(mplc.cnames[s])  # a tuple
                        #     print(rgb_tuple)
                        #     rgb = [x for x in rgb_tuple]
                        #     print(rgb)
    
                    #    if (child.attrib['id'] == 'cells'):
                    #      print('-------- found cells!!')
                    #      tissue_child = child
    
                    #### num_cells becomes total number of cells per frame/sample
                    num_cells += 1
                print(fname, ':  num_cells= ', num_cells)

            return d, d_attributes, title_str, plot_x_extend, plot_y_extend

    def create_cell_layer_from_SVG(self, cell_positions: dict, cell_attributes: dict):
        d = cell_positions
        d_attributes = cell_attributes

        ####################################################################################################################
        ####################################        Plot cell tracks and other options              ########################
        ####################################################################################################################

        # ax.set_xticks([])
        # ax.set_yticks([]);
        # ax.set_xlim(0, 8); ax.set_ylim(0, 8)

        # print 'dir(fig)=',dir(fig)
        # fig.set_figwidth(8)
        # fig.set_figheight(8)

        count = 0

        # weighting = np.linspace(0.0001, 3.5, num=number_of_samples)
        #
        # weighting = np.log10(weighting)

        ##### Extract and plot position data for each cell found
        for key in d.keys():
            if (len(d[key].shape) == 2):
                x = d[key][:, 0]
                y = d[key][:, 1]

                # plt.plot(x, y,'-')  # plot doesn't seem to allow weighting or size variation at all in teh connections ...  # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.arrow.html or https://stackoverflow.com/questions/7519467/line-plot-with-arrows-in-matplotlib
                # plt.scatter(x, y, s = weighting) - scatter allows weighting but doens't connect ...
                # plt.scatter(x, y, s=weighting) # could try a non-linear weighting ...

                #### Plot cell track as a directed, weighted (by length) path
                self.ax.quiver(x[:-1], y[:-1], x[1:] - x[:-1], y[1:] - y[:-1], scale_units='xy', angles='xy', scale=1,
                               minlength=0.001, headwidth=1.5, headlength=4)

                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                # self.ax.scatter(x[-1], y[-1], s=85.0, c=d_attributes[key], alpha=0.7)

                # Add cells layer
                # for i, ct in enumerate(types):
                #     plot_df = cell_df[cell_df['cell_type'] == ct]
                #     for j in plot_df.index:
                circ = Circle((x[-1], y[-1]),
                                      radius=8.41271, color=d_attributes[key], alpha=0.7)
                        # for a blue circle with a black edge
                        # circ = Circle((plot_df.loc[j, 'position_x'], plot_df.loc[j, 'position_y']),
                        #                radius=plot_df.loc[j, 'radius'], alpha=0.7, edgecolor='black')
                self.ax.add_artist(circ)

            #### used if history lenght is set to 0 and if in first frame of sequnece (there is no history)
            elif (len(d[key].shape) == 1):
                x = d[key][0]
                y = d[key][1]
                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                circ = Circle((x, y),
                              radius=8.41271, color=d_attributes[key], alpha=0.7)
                self.ax.add_artist(circ)
                # self.ax.scatter(x, y, s=85.0, c=d_attributes[key], alpha=0.7)
                # plt.scatter(x, y, s=3.5, c=)

            else:
                print(key, " has no x,y points")

    def create_figure_from_SVG (self, cell_positions: dict, cell_attributes: dict):
        d = cell_positions
        d_attributes = cell_attributes
    
    
        ####################################################################################################################
        ####################################        Plot cell tracks and other options              ########################
        ####################################################################################################################
    
        # ax.set_xticks([])
        # ax.set_yticks([]);
        # ax.set_xlim(0, 8); ax.set_ylim(0, 8)
    
        # print 'dir(fig)=',dir(fig)
        # fig.set_figwidth(8)
        # fig.set_figheight(8)
    
        count = 0
    
        # weighting = np.linspace(0.0001, 3.5, num=number_of_samples)
        #
        # weighting = np.log10(weighting)
    
        ##### Extract and plot position data for each cell found
        for key in d.keys():
            if (len(d[key].shape) == 2):
                x = d[key][:, 0]
                y = d[key][:, 1]
    
                # plt.plot(x, y,'-')  # plot doesn't seem to allow weighting or size variation at all in teh connections ...  # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.arrow.html or https://stackoverflow.com/questions/7519467/line-plot-with-arrows-in-matplotlib
                # plt.scatter(x, y, s = weighting) - scatter allows weighting but doens't connect ...
                # plt.scatter(x, y, s=weighting) # could try a non-linear weighting ...
    
                #### Plot cell track as a directed, weighted (by length) path
                self.ax.quiver(x[:-1], y[:-1], x[1:] - x[:-1], y[1:] - y[:-1], scale_units='xy', angles='xy', scale=1,
                           minlength=0.001, headwidth=1.5, headlength=4)
    
                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                self.ax.scatter(x[-1], y[-1], s=85.0, c=d_attributes[key], alpha=0.7)
    
            #### used if history lenght is set to 0 and if in first frame of sequnece (there is no history)
            elif (len(d[key].shape) == 1):
                x = d[key][0]
                y = d[key][1]
                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                self.ax.scatter(x, y, s=85.0, c=d_attributes[key], alpha=0.7)
                # plt.scatter(x, y, s=3.5, c=)
    
            else:
                print(key, " has no x,y points")

    def produce_movie(self, data_path: str=  '.',  save_path: str= '', start_file_index: int = 0, sample_step_interval: int = 1, 
                             end_file_index: int=120, trail_length: int=10, movie_options: dict=None, image_options: dict=None):
        if movie_options is None:
            movie_options =  {'INCLUDE_ALL_SVGs': True,
                            'INCLUDE_FULL_HISTORY': True
                            }

        if image_options is None:
            image_options = {"produce_for_movie" : True,
                        "show_plot": False}
            # movie_options['INCLUDE_ALL_SVGs'] = True
            # movie_options['INCLUDE_FULL_HISTORY'] = True

            #### Get list of all file names in directory

        # data_path: str, save_path: str, save_name: str, start_file_index: int, end_file_index: int,
                # trail_length: int, INCLUDE_ALL_SVGs: bool, INCLUDE_FULL_HISTORY: bool)

        # def generic_plotter(self, starting_index: int = 0, sample_step_interval: int = 1, number_of_samples: int = 120,
        #             file_name: str = None, input_path: str= '.', output_path: str= '', naming_index: int=0, options=None):
        
        files = os.listdir(data_path)

        list_of_svgs = []

        #### examine all file names in directory and add ones, via string matching, as needed to list of names of files of interest
        for i in range(len(files)):
            if not re.search('snapshot(.*)\.svg', files[i]):
                continue

            # I feel like a dictionary could be used here, but I really need some ordering. A dict might be faster, but I don't
            # expect huge file lists. So I will just sort as I know how to do that ...

            list_of_svgs.append(files[i])

        #### Sort file name list
        list_of_svgs.sort()

        truncated_list_of_svgs = []

        #### Reduce file list to times of interst only
        for i in range(len(list_of_svgs)):

            if i < start_file_index:
                continue

            if i >= end_file_index:
                continue

            truncated_list_of_svgs.append(list_of_svgs[i])

        # print(list_of_svgs)
        print(truncated_list_of_svgs)

        if movie_options['INCLUDE_ALL_SVGs'] :
            print('Including all SVGs')
            truncated_list_of_svgs = list_of_svgs

        max_number_of_samples = trail_length

        if movie_options['INCLUDE_FULL_HISTORY']:
            print('Including full positional history of cells')
            max_number_of_samples = len(truncated_list_of_svgs)

        print('Processing {} SVGs'.format(len(truncated_list_of_svgs)))

        # Also, as written it isn't very flexible
        # would certainly be ideal to not call plot_cell_tracks every time, but instead store what is available. Could add a function that just
        # extracts the data from one SVG then appends it to exsisting data structure. could read all the desired data into Pandas DF
        # then write out images. Etc. But as is, this is definitely reading the SVGs much to frequently.

        for i in range(len(truncated_list_of_svgs)):
            j = i + 1 # this offsets the index so that we don't report that 0 samples have been taken, while stil producing an image.
            starting_index = j - max_number_of_samples

            #### Goes with "trail closing" block - not currently being used.
            projected_upper_sample_index = max_number_of_samples + starting_index
            max_samples_left = len(truncated_list_of_svgs) - j

            # def generic_plotter(self, starting_index: int = 0, sample_step_interval: int = 1, number_of_samples: int = 120,
    #             file_name: str = None, input_path: str= '.', output_path: str= '', naming_index: int=0, options=None):

            if i >= max_number_of_samples:
                self.generic_plotter(starting_index, 1, max_number_of_samples, naming_index=i, options=image_options)
                # print('middle')

            #### If one wanted to make the trails collapse into the last available location of the cell you would use something
            #### like this elif block
            # elif projected_upper_sample_index > len(list_of_svgs)-1:
            #     plot_cell_tracks(starting_index, 1, max_samples_left, True, True, i)
            #     print(max_samples_left)
            #     print('late')
            else:
                self.generic_plotter(0, 1, j, naming_index=i, options=image_options)
                # print('early')

        #### Total frames to include in moview
        number_frames = end_file_index - start_file_index
        
        if movie_options['INCLUDE_ALL_SVGs']:
            number_frames = len(list_of_svgs)
            start_file_index = 0

        # string_of_interest = 'ffmpeg -start_number ' + str(
        #     start_file_index) + ' -y -framerate 12 -i ' + save_path + 'output%08d.png' + ' -frames:v ' + str(
        #     number_frames) + ' -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"'
        # print(string_of_interest)
        os.system(
            'ffmpeg -start_number ' + str(
                start_file_index) + ' -y -framerate 12 -i ' + save_path + 'output%08d.png' + ' -frames:v ' + str(
                number_frames) + ' -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"')

        # https://superuser.com/questions/666860/clarification-for-ffmpeg-input-option-with-image-files-as-input
        # https://superuser.com/questions/734976/ffmpeg-limit-number-of-images-converted-to-video


    def general_image_plotter (filename: str=None, folder: str='.', output_folder='', cell_df: dict=None, cell_positions_from_SVG: dict=None, cell_attributes_from_SVG: dict=None, chemical_mesh: dict=None, ECM_mesh: dict=None, options=None):
        if options is None:
            options = {"output_plot": True,
                       "show_plot": True,
                       "produce_for_panel": False
                       }
        output_plot = options['output_plot']
        show_plot = options['show_plot']
        produce_for_panel = options['produce_for_panel']
    
        # Testing this
    
        d = cell_positions_from_SVG
        d_attributes = cell_attributes_from_SVG
        
        ####################################################################################################################
        ####################################        Plot cell tracks and other options              ########################
        ####################################################################################################################
    
        fig = plt.figure(figsize=(7, 7))
        ax = fig.gca()
        ax.set_aspect("equal")
        # ax.set_xticks([])
        # ax.set_yticks([]);
        # ax.set_xlim(0, 8); ax.set_ylim(0, 8)
    
        # print 'dir(fig)=',dir(fig)
        # fig.set_figwidth(8)
        # fig.set_figheight(8)
    
        count = 0
    
        # weighting = np.linspace(0.0001, 3.5, num=number_of_samples)
        #
        # weighting = np.log10(weighting)
    
        ##### Extract and plot position data for each cell found
        for key in d.keys():
            if (len(d[key].shape) == 2):
                x = d[key][:, 0]
                y = d[key][:, 1]
    
                # plt.plot(x, y,'-')  # plot doesn't seem to allow weighting or size variation at all in teh connections ...  # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.arrow.html or https://stackoverflow.com/questions/7519467/line-plot-with-arrows-in-matplotlib
                # plt.scatter(x, y, s = weighting) - scatter allows weighting but doens't connect ...
                # plt.scatter(x, y, s=weighting) # could try a non-linear weighting ...
    
                #### Plot cell track as a directed, weighted (by length) path
                plt.quiver(x[:-1], y[:-1], x[1:] - x[:-1], y[1:] - y[:-1], scale_units='xy', angles='xy', scale=1,
                           minlength=0.001, headwidth=1.5, headlength=4)
    
                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                plt.scatter(x[-1], y[-1], s=85.0, c=d_attributes[key], alpha=0.7)
    
            #### used if history lenght is set to 0 and if in first frame of sequnece (there is no history)
            elif (len(d[key].shape) == 1):
                x = d[key][0]
                y = d[key][1]
                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                plt.scatter(x, y, s=85.0, c=d_attributes[key], alpha=0.7)
                # plt.scatter(x, y, s=3.5, c=)
    
            else:
                print(key, " has no x,y points")
    
        #### Build plot frame, titles, and save data
        plt.ylim(-1000 / 2, 1000 / 2)
        plt.xlim(-1000 / 2, 1000 / 2)
    
        if produce_for_panel == False:
            title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str
            # %"Starting at frame {}, sample interval of {} for {} total samples".format(number_of_samples, sample_step_interval, number_of_samples)
            plt.title(title_str)
        else:
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            ax.set_xlabel('microns', fontsize=20)
            ax.set_ylabel('microns', fontsize=20)
            fig.tight_layout()
        # could change to the custom in the movie output or some other more better output if desired.
        output_folder = ''
        if file_name is None:
            file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
    
        # Produce plot following the available options.
    
        if output_plot is True:
            plt.savefig(output_folder + file_name + '.png', dpi=256)
        if show_plot is True:
            plt.show()
        # plt.close()
        return fig
    
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    
    def plot_cells_and_uE_for_movie (starting_index: int, sample_step_interval: int, number_of_samples: int, naming_index: int, options=None ):
        if options is None:
            options= {"output_plot": True,
                      "show_plot": False,
                      "produce_for_panel": True
                      }
    
        ### Now add in place to call the regular one, but with a string ...
        file_name = 'output' + f'{naming_index:08}'
        print(options)
        plot_cell_tracks_from_svg (starting_index, sample_step_interval, number_of_samples, file_name, options)
        # or a dictionary - and then I just modify the dictionary for the options or even have several of them
    
    
        # 0
        # 1
        # 417
        # True
        # True
        # True
    
    # test_of_args_and_kwargs(0, 0, 0, output_plot='False')
    
    def plot_cell_tracks_from_svg(starting_index: int, sample_step_interval: int, number_of_samples: int, file_name: str=None, options=None ):
        """
        Produces savable image of cell positional history, plotted as arrows (quiver plot) with final cell positions plotted as a cirle.
        Slight modification of the function in cell_track_plotter. The modification allows for tracking the index of a series
        of inputs such that outputs of this function can be appropriate indexed and compiled into a movie.
    
        sample_step_interval * number_of_samples - starting_index yields the trail length in time steps. number_of_samples provides
        the number of intervals plotted per image.
    
        Parameters
        ----------
        starting_index :
            Integer index of the PhysiCell SVG output to begin trackign at
        sample_step_interval :
            Interval (number of time steps (SVGs)) to sample at. A value of 2 would add a tracking point for every other SVG
        number_of_samples :
            Number of SVGs to process (total)/Length of cell positional history. Number of samples * sample size step interval provides the index of the final SVG to process
        output_plot :
            Save plot flag (required to produce a movie from resulting images)
        show_plot :
            Show plot flag (for processing many images in a loop, this should likely be set to false. Images have to be closed manually)
        produce_for_panel :
                Flag - calls tight_layout, increases axes font sizes, and plots without title. For using in panels of images where there will be captions.
        Returns
        -------
        Null :
            Produces a png image from the input PhysiCell SVGs.
    
        """
    
        if options is None:
            options= {"output_plot": True,
                      "show_plot": True,
                      "produce_for_panel": False
                      }
        output_plot = options['output_plot']
        show_plot = options['show_plot']
        produce_for_panel = options['produce_for_panel']
    
        d={}   # dictionary to hold all (x,y) positions of cells
        d_attributes = {}   #dictionary to hold other attributes, like color (a data frame might be nice here in the long run ... ) \
                            # currently only being read once only as cell dictionary is populated - so only use for static values!
    
        """ 
        --- for example ---
        In [141]: d['cell1599'][0:3]
        Out[141]: 
        array([[ 4900.  ,  4900.  ],
               [ 4934.17,  4487.91],
               [ 4960.75,  4148.02]])
        """
    
        ####################################################################################################################
        ####################################            Generate list of file indices to load       ########################
        ####################################################################################################################
    
        endpoint = starting_index + sample_step_interval*number_of_samples
        file_indices = np.linspace(starting_index, endpoint, num=number_of_samples, endpoint=False)
        print(file_indices)
    
        maxCount = starting_index
    
        ####### Uncomment for statement below to generate a random list of file names versus the prespecifed list. ########
        ####### Leaving for historical record. If used, the inputs would need to be a single integer,       ########
        ####### versus the three integers required to generate the prespecified list. Also, remove the other for statement. ########
        # count = 0
        #
        # for fname in glob.glob('snapshot*.svg'):
        #     print(fname)
        # # for fname in['snapshot00000000.svg', 'snapshot00000001.svg']:
        # # for fname in['snapshot00000000.svg']:
        # #   print(fname)
        #   count += 1
        #   if count > maxCount:
        #     break
    
    
        ####################################################################################################################
        ####################################        Main loading and processing loop                ########################
        ####################################################################################################################
    
        for file_index in file_indices:
            fname = "%0.8d" % file_index
            fname = 'snapshot' + fname + '.svg'# https://realpython.com/python-f-strings/
            print(fname)
    
            ##### Parse XML tree into a dictionary called 'tree" and get root
            # print('\n---- ' + fname + ':')
            tree=ET.parse(fname)
    
            # print('--- root.tag, root.attrib ---')
            root=tree.getroot()
            #  print('--- root.tag ---')
            #  print(root.tag)
            #  print('--- root.attrib ---')
            #  print(root.attrib)
            #  print('--- child.tag, child.attrib ---')
    
            numChildren = 0
    
            ### Find branches coming from root - tissue parents
            for child in root:
                # print(child.tag, child.attrib)
                #    print('attrib=',child.attrib)
                #  if (child.attrib['id'] == 'tissue'):
    
                if child.text and "Current time" in child.text:
                    svals = child.text.split()
                    title_str = "Current time: " + svals[2] + "d, " + svals[4] + "h, " + svals[
                        7] + "m"
    
                if 'width' in child.attrib.keys():
                    #### Assumes a 70 length unit offsite inthe the Y dimension of the SVG!!!!!!
                    plot_x_extend = float(child.attrib['width'])
                    plot_y_extend = float(child.attrib['height'])
    
                    #### Remove the padding placed into the SVG to determine the true y extend
                    plot_y_extend = plot_y_extend-70
    
                    #### Find the coordinate transform amounts
                    y_coordinate_transform = plot_y_extend/2
                    x_coordinate_transform = plot_x_extend/2
    
                ##### Find the tissue tag and make it child
                if 'id' in child.attrib.keys():
                #      print('-------- found tissue!!')
                    tissue_parent = child
                    break
    
            #  print('------ search tissue')
    
            ### find the branch with the cells "id=cells" among all the branches in the XML root
            for child in tissue_parent:
            #    print('attrib=',child.attrib)
                if (child.attrib['id'] == 'cells'):
                #      print('-------- found cells!!')
                    cells_parent = child
                    break
                numChildren += 1
    
            ### Search within the cells branch for all indiviual cells. Get their locations
            num_cells = 0
            #  print('------ search cells')
            for child in cells_parent:
            #    print(child.tag, child.attrib)
            #    print('attrib=',child.attrib)
    
                # Find the locations of the cells within the cell tags
                for circle in child:
                #      print('  --- cx,cy=',circle.attrib['cx'],circle.attrib['cy'])
                    xval = float(circle.attrib['cx'])
    
                    # should we test for bogus x,y locations??
                    if (math.fabs(xval) > 10000.):
                        print("xval=",xval)
                        break
                    yval = float(circle.attrib['cy']) #- y_coordinate_transform
                    if (math.fabs(yval) > 10000.):
                        print("yval=",yval)
                        break
    
                    # Pull out the cell's location. If ID not already in stack to track, put in new cell in dictionary while applying coordinate transform.
                    if (child.attrib['id'] in d.keys()):
                        d[child.attrib['id']] = np.vstack((d[child.attrib['id']], [ float(circle.attrib['cx'])-x_coordinate_transform, float(circle.attrib['cy'])-y_coordinate_transform ]))
                    #### Comment out this else to produce single cell tracks
                    else:
                        d[child.attrib['id']] = np.array( [ float(circle.attrib['cx'])-x_coordinate_transform, float(circle.attrib['cy'])-y_coordinate_transform])
                        d_attributes[child.attrib['id']] = circle.attrib['fill']
    
                    ###### Uncomment this elif and else to produce single cell tracks
                    # elif (child.attrib['id'] == 'cell24'):
                    #     d[child.attrib['id']] = np.array( [ float(circle.attrib['cx'])-x_coordinate_transform, float(circle.attrib['cy'])-y_coordinate_transform])
                    #     d_attributes[child.attrib['id']] = circle.attrib['fill']
                    # else:
                    #     break
    
                    ##### This 'break' statement is required to skip the nucleus circle. There are two circle attributes. \
                    ##### If both nuclear and cell boundary attributes are needed, this break NEEDS REMOVED!!!!
                    break
    
                    ### Code to translate string based coloring to rgb coloring. Use as needed.
                      # s = circle.attrib['fill']
                      # print("s=",s)
                      # print("type(s)=",type(s))
                      # if (s[0:3] == "rgb"):  # if an rgb string, e.g. "rgb(175,175,80)"
                      #     #  circle.attrib={'cx': '1085.59','cy': '1225.24','fill': 'rgb(159,159,96)','r': '6.67717','stroke': 'rgb(159,159,96)','stroke-width': '0.5'}
                      #     rgb = list(map(int, s[4:-1].split(",")))
                      #     rgb[:] = [x / 255. for x in rgb]
                      # else:  # otherwise, must be a color name
                      #     rgb_tuple = mplc.to_rgb(mplc.cnames[s])  # a tuple
                      #     print(rgb_tuple)
                      #     rgb = [x for x in rgb_tuple]
                      #     print(rgb)
    
            #    if (child.attrib['id'] == 'cells'):
            #      print('-------- found cells!!')
            #      tissue_child = child
    
                #### num_cells becomes total number of cells per frame/sample
                num_cells += 1
            print(fname,':  num_cells= ',num_cells)
    
        ####################################################################################################################
        ####################################        Plot cell tracks and other options              ########################
        ####################################################################################################################
    
        fig = plt.figure(figsize=(7,7))
        ax = fig.gca()
        ax.set_aspect("equal")
        #ax.set_xticks([])
        #ax.set_yticks([]);
        #ax.set_xlim(0, 8); ax.set_ylim(0, 8)
    
        #print 'dir(fig)=',dir(fig)
        #fig.set_figwidth(8)
        #fig.set_figheight(8)
    
        count = 0
    
        # weighting = np.linspace(0.0001, 3.5, num=number_of_samples)
        #
        # weighting = np.log10(weighting)
    
        ##### Extract and plot position data for each cell found
        for key in d.keys():
            if (len(d[key].shape) == 2):
                x = d[key][:,0]
                y = d[key][:,1]
    
                # plt.plot(x, y,'-')  # plot doesn't seem to allow weighting or size variation at all in teh connections ...  # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.arrow.html or https://stackoverflow.com/questions/7519467/line-plot-with-arrows-in-matplotlib
                # plt.scatter(x, y, s = weighting) - scatter allows weighting but doens't connect ...
                # plt.scatter(x, y, s=weighting) # could try a non-linear weighting ...
    
                #### Plot cell track as a directed, weighted (by length) path
                plt.quiver(x[:-1], y[:-1], x[1:] - x[:-1], y[1:] - y[:-1], scale_units='xy', angles='xy', scale=1, minlength = 0.001, headwidth=1.5, headlength=4)
    
                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                plt.scatter(x[-1], y[-1], s=85.0, c=d_attributes[key], alpha=0.7)
    
            #### used if history lenght is set to 0 and if in first frame of sequnece (there is no history)
            elif (len(d[key].shape) == 1):
                x = d[key][0]
                y = d[key][1]
                #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
                plt.scatter(x, y, s=85.0, c=d_attributes[key], alpha=0.7)
                # plt.scatter(x, y, s=3.5, c=)
    
            else:
                print(key, " has no x,y points")
    
        #### Build plot frame, titles, and save data
        plt.ylim(-plot_y_extend/2, plot_y_extend/2)
        plt.xlim(-plot_x_extend/2, plot_x_extend/2)
    
        if produce_for_panel == False:
            title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str
            # %"Starting at frame {}, sample interval of {} for {} total samples".format(number_of_samples, sample_step_interval, number_of_samples)
            plt.title(title_str)
        else:
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            ax.set_xlabel('microns', fontsize=20)
            ax.set_ylabel('microns', fontsize=20)
            fig.tight_layout()
        # could change to the custom in the movie output or some other more better output if desired.
        output_folder = ''
        if file_name is None:
            file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
    
        # Produce plot following the available options.
    
        if output_plot is True:
            plt.savefig(output_folder + file_name + '.png', dpi=256)
        if show_plot is True:
            plt.show()
        # plt.close()
    
    def cell_history_movie_generator_from_SVG(data_path: str, save_path: str, save_name: str, start_file_index: int, end_file_index: int,
                            trail_length: int, INCLUDE_ALL_SVGs: bool, INCLUDE_FULL_HISTORY: bool):
    
        """
        Generates the list of files in data_path, finds the relevant SVGs, makes plots from them, then outputs an
        ffmpeg generated movie to save_path, naming the movie save_name.
    
        This function requires ffmpeg be installed at the command line.
    
    
        :param data_path: Path to directory containing data
        :param save_path: Path to save generated image(s) and movie to
        :param save_name: Save name for movie
        :param start_file_index: For the plotting call - Integer index of the PhysiCell SVG output to begin tracking at
        :param end_file_index: For the plotting call - Integer index of last PhysiCell SVG output to include in movie
        :param trail_length: For the plotting call - Length (in output steps) of cell positional history to include in movie
        :param INCLUDE_ALL_SVGs: If true, all findable PhysiCell SVGs are processed and included in movie
        :param INCLUDE_FULL_HISTORY: If true, the entire available cell history is included, regardless of the value of trail length.
        :return: Null. Produces a series of images from PhysiCell SVGs and movie from said images.
        """
    
        #### Flags (for cell track plotter calls)
    
        output_plot = True
        show_plot = False
        produce_for_panel = False
    
        #### Get list of all file names in directory
        files = os.listdir(data_path)
    
        list_of_svgs = []
    
        #### examine all file names in directory and add ones, via string matching, as needed to list of names of files of interest
        for i in range(len(files)):
            if not re.search('snapshot(.*)\.svg', files[i]):
                continue
    
            # I feel like a dictionary could be used here, but I really need some ordering. A dict might be faster, but I don't
            # expect huge file lists. So I will just sort as I know how to do that ...
    
            list_of_svgs.append(files[i])
    
        #### Sort file name list
        list_of_svgs.sort()
    
        truncated_list_of_svgs = []
    
        #### Reduce file list to times of interst only
        for i in range(len(list_of_svgs)):
    
            if i < start_file_index:
                continue
    
            if i >= end_file_index:
                continue
    
            truncated_list_of_svgs.append(list_of_svgs[i])
    
        # print(list_of_svgs)
        print(truncated_list_of_svgs)
    
        if INCLUDE_ALL_SVGs:
            print('Including all SVGs')
            truncated_list_of_svgs = list_of_svgs
    
        max_number_of_samples = trail_length
    
        if INCLUDE_FULL_HISTORY:
            print('Including full positional history of cells')
            max_number_of_samples = len(truncated_list_of_svgs)
    
        print('Processing {} SVGs'.format(len(truncated_list_of_svgs)))
    
        # Also, as written it isn't very flexible
        # would certainly be ideal to not call plot_cell_tracks every time, but instead store what is available. Could add a function that just
        # extracts the data from one SVG then appends it to exsisting data structure. could read all the desired data into Pandas DF
        # then write out images. Etc. But as is, this is definitely reading the SVGs much to frequently.
    
        for i in range(len(truncated_list_of_svgs)):
            j = i + 1  # this offsets the index so that we don't report that 0 samples have been taken, while stil producing an image.
            starting_index = j - max_number_of_samples
    
            #### Goes with "trail closing" block - not currently being used.
            projected_upper_sample_index = max_number_of_samples + starting_index
            max_samples_left = len(truncated_list_of_svgs) - j
    
            if i >= max_number_of_samples:
                plot_cell_tracks_for_movie(starting_index, 1, max_number_of_samples, output_plot, show_plot, i,
                                           produce_for_panel)
                # print('middle')
    
            #### If one wanted to make the trails collapse into the last available location of the cell you would use something
            #### like this elif block
            # elif projected_upper_sample_index > len(list_of_svgs)-1:
            #     plot_cell_tracks(starting_index, 1, max_samples_left, True, True, i)
            #     print(max_samples_left)
            #     print('late')
            else:
                plot_cell_tracks_for_movie(0, 1, j, output_plot, show_plot, i, produce_for_panel)
                # print('early')
    
        #### Total frames to include in moview
        number_frames = end_file_index - start_file_index
    
        if INCLUDE_ALL_SVGs:
            number_frames = len(list_of_svgs)
            start_file_index = 0
    
        # string_of_interest = 'ffmpeg -start_number ' + str(
        #     start_file_index) + ' -y -framerate 12 -i ' + save_path + 'output%08d.png' + ' -frames:v ' + str(
        #     number_frames) + ' -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"'
        # print(string_of_interest)
        os.system(
            'ffmpeg -start_number ' + str(
                start_file_index) + ' -y -framerate 12 -i ' + save_path + 'output%08d.png' + ' -frames:v ' + str(
                number_frames) + ' -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"')
    
        # https://superuser.com/questions/666860/clarification-for-ffmpeg-input-option-with-image-files-as-input
        # https://superuser.com/questions/734976/ffmpeg-limit-number-of-images-converted-to-video


