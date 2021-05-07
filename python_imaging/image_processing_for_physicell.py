import sys
import xml.etree.ElementTree as ET
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import math
import distutils.util
import pyMCDS_ECM

# pseudo code:

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

class PhysiCellPlotter():
    
    def __init__(self, parent = None):
        self.figsize_width_svg = 7.0
        self.figsize_height_svg = 7.0
        self.title = "title"
        self.fig, self.ax = plt.subplots(figsize=(self.figsize_width_svg, self.figsize_height_svg))

    def generic_plotter(self, starting_index: int = 0, sample_step_interval: int = 1, number_of_samples: int = 120,
                        file_name: str = None, options=None):
        #### Needs output and input folders!!!!!!!!

        if options is None:
            options = {"output_plot": True,
                       "show_plot": True,
                       "produce_for_panel": False
                       }
        output_plot = options['output_plot']
        show_plot = options['show_plot']
        produce_for_panel = options['produce_for_panel']

        cell_positions, cell_attributes, title_str, plot_x_extend, plot_y_extend = self.load_cell_positions_from_SVG(
            starting_index, sample_step_interval, number_of_samples)
        # print(type(cell_attributes))
        # print(type(cell_positions))
        # print(title_str, plot_x_extend)
        # fig = plt.figure(figsize=(7, 7))
        # print(fig)
        # self.fig.show()
        self.create_figure_from_SVG(cell_positions, cell_attributes)
        # fig = plt.gcf()
        # plt.show(figsize=(7, 7))

        endpoint = starting_index + sample_step_interval * number_of_samples
        title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str
        if file_name is None:
            file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)

        self.plot_figure(title_str, plot_x_extend, plot_y_extend, file_name, options)

        # self.ax.set_aspect("equal")
        # # endpoint = starting_index + sample_step_interval*number_of_samples
        # #### Build plot frame, titles, and save data
        # self.ax.set_ylim(-plot_y_extend / 2, plot_y_extend / 2)
        # self.ax.set_xlim(-plot_x_extend / 2, plot_x_extend / 2)

        # if produce_for_panel == False:
        #     # title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str
        #     # %"Starting at frame {}, sample interval of {} for {} total samples".format(number_of_samples, sample_step_interval, number_of_samples)
        #     self.ax.set_title(title_str)
        # else:
        #     self.ax.xticks(fontsize=20)
        #     self.ax.yticks(fontsize=20)
        #     self.ax.set_xlabel('microns', fontsize=20)
        #     self.ax.set_ylabel('microns', fontsize=20)
        #     self.fig.tight_layout()
        # # could change to the custom in the movie output or some other more better output if desired.
        # output_folder = ''
        # # if file_name is None:
        # #     file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
        #
        # # Produce plot following the available options.
        #
        # if output_plot is True:
        #     plt.savefig(output_folder + file_name + '.png', dpi=256)
        # if show_plot is True:
        #     plt.show()

        # fig.figure(figsize=(7, 7))
        # ax = fig.gca()
        # ax.set_aspect("equal")
        # endpoint = starting_index + sample_step_interval*number_of_samples
        # #### Build plot frame, titles, and save data
        # plt.ylim(-plot_y_extend/2, plot_y_extend/2)
        # plt.xlim(-plot_x_extend/2, plot_x_extend/2)
        #
        # if produce_for_panel == False:
        #     title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str
        #     # %"Starting at frame {}, sample interval of {} for {} total samples".format(number_of_samples, sample_step_interval, number_of_samples)
        #     plt.title(title_str)
        # else:
        #     plt.xticks(fontsize=20)
        #     plt.yticks(fontsize=20)
        #     ax.set_xlabel('microns', fontsize=20)
        #     ax.set_ylabel('microns', fontsize=20)
        #     fig.tight_layout()
        # # could change to the custom in the movie output or some other more better output if desired.
        # output_folder = ''
        # if file_name is None:
        #     file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
        #
        # # Produce plot following the available options.
        #
        # if output_plot is True:
        #     plt.savefig(output_folder + file_name + '.png', dpi=256)
        # if show_plot is True:
        #     plt.show()

        # return d, d_attributes, title_str, plot_x_extend, plot_y_extend
    
    def plot_figure(self, title_str, plot_x_extend, plot_y_extend, file_name, options: dict=None): ###### This should probably have to have options??????? Why though???
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
            self.ax.xticks(fontsize=20)
            self.ax.yticks(fontsize=20)
            self.ax.set_xlabel('microns', fontsize=20)
            self.ax.set_ylabel('microns', fontsize=20)
            self.fig.tight_layout()
        # could change to the custom in the movie output or some other more better output if desired.
        output_folder = ''
        # if file_name is None:
        #     file_name = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
    
        # Produce plot following the available options.
    
        if output_plot is True:
            plt.savefig(output_folder + file_name + '.png', dpi=256)
        if show_plot is True:
            plt.show()
    
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
                        #### Assumes a 70 length unit offsite inthe the Y dimension of the SVG!!!!!!
                        plot_x_extend = float(child.attrib['width'])
                        plot_y_extend = float(child.attrib['height'])
    
                        #### Remove the padding placed into the SVG to determine the true y extend
                        plot_y_extend = plot_y_extend - 70
    
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


