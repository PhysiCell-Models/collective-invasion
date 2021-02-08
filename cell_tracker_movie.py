#
# cell_tracker.py - plot 2-D cell tracks associated with PhysiCell .svg files
#
# Usage:
#  Takes 0, 1, or 7 arguments. See below line 239 in "if __name__ == '__main__':" for usage.
#
# Dependencies include matplotlib and numpy. We recommend installing the Anaconda Python3 distribution.
#
# Examples (run from directory containing the .svg files):
#  See below line 239 in "if __name__ == '__main__':"
#
# Author: function plot_cell_tracks_for_movie - Randy Heiland, modified by John Metzcar (see cell_track_plotter.py and cell_tracks.py as well for original functions)
#         This script cell_tracker_movie.py - John Metzcar (Twitter - @jmetzcar). See also anim_svg_opac.py in PhysiCell tools for coloring functionality

import sys
import xml.etree.ElementTree as ET
import numpy as np
import glob
import matplotlib.pyplot as plt
import math, os, sys, re
import distutils.util

def plot_cell_tracks_for_movie(starting_index: int, sample_step_interval: int, number_of_samples: int, output_plot: bool,
                     show_plot: bool, naming_index: int, produce_for_panel: bool):
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
        Interval (number of time steps (SVGs)) to sample at. A value of 2 would add a tracking point for every other SVG. For this special
        function, it is currently (01.27.21) assumed that will be 1.
    number_of_samples :
        Number of SVGs to process (total)/Length of cell positional history. Number of samples * sample size step interval provides the index of the final SVG to process
    output_plot :
        Save plot flag (required to produce a movie from resulting images)
    show_plot :
        Show plot flag (for processing many images in a loop, this should likely be set to false. Images have to be closed manually)
    naming_index :
        Unique to this function. Index used in naming output file of plot_cell_tracks function - filename = output + naming_index.png and leading zeros as needed.
    produce_for_panel :
        Flag - calls tight_layout, increases axes font sizes, and plots without title. For using in panels of images where there will be captions.
    Returns
    -------
    Null :
        Produces a png image from the input PhysiCell SVGs.
    """

    #### Flags

    output_plot = output_plot
    show_plot = show_plot
    naming_index = naming_index
    produce_for_panel = produce_for_panel

    d = {}  # dictionary to hold all (x,y) positions of cells
    d_attributes = {}   #dictionary to hold other attributes, like color (a data frame might be nice here in the long run ... )\
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

    ####################################################################################################################
    ####################################        Main loading and processing loop                ########################
    ####################################################################################################################

    for file_index in file_indices:
        fname = "%0.8d" % file_index
        fname = 'snapshot' + fname + '.svg'  # https://realpython.com/python-f-strings/
        # print(fname)

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
            ##### Find the tissue tag and make it child

            if child.text and "Current time" in child.text:
                svals = child.text.split()
                title_str = "Current time: " + svals[2] + "d, " + svals[4] + "h, " + svals[
                    7] + "m"

            if ('width' in child.attrib.keys()):
                #### Assumes a 70 length unit offsite inthe the Y dimension of the SVG!!!!!!
                plot_x_extend = float(child.attrib['width'])
                plot_y_extend = float(child.attrib['height'])

                #### Remove the padding placed into the SVG to determine the true y extend
                plot_y_extend = plot_y_extend-70

                #### Find the coordinate transform amounts
                y_coordinate_transform = plot_y_extend/2
                x_coordinate_transform = plot_x_extend/2

            if ('id' in child.attrib.keys()):
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
                yval = float(circle.attrib['cy'])
                if (math.fabs(yval) > 10000.):
                    print("xval=", xval)
                    break

                # Pull out the cell's location. If ID not already in stack to track, put in new cell in dictionary
                if (child.attrib['id'] in d.keys()):
                    d[child.attrib['id']] = np.vstack((d[child.attrib['id']],
                                                       [float(circle.attrib['cx']) - x_coordinate_transform,
                                                        float(circle.attrib['cy']) - y_coordinate_transform]))
                else:
                    d[child.attrib['id']] = np.array([float(circle.attrib['cx']) - x_coordinate_transform,
                                                      float(circle.attrib['cy']) - y_coordinate_transform])
                    d_attributes[child.attrib['id']] = circle.attrib['fill']
                ##### This 'break' statement is required to skip the nucleus circle. There are two circle attributes. \
                ##### If both nuclear and cell boundary attributes are needed, this break NEEDS REMOVED!!!!
                break

            #    if (child.attrib['id'] == 'cells'):
            #      print('-------- found cells!!')
            #      tissue_child = child

            #### num_cells becomes total number of cells per frame/sample
            num_cells += 1
        print(fname, ':  num_cells= ', num_cells)

    ####################################################################################################################
    ####################################        Plot cell tracks and other options              ########################
    ####################################################################################################################

    fig = plt.figure(figsize=(8, 8))
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
            plt.scatter(x[-1], y[-1], s = 5.0, c = d_attributes[key])

        #### used if history lenght is set to 0 and if in first frame of sequnece (there is no history)
        elif (len(d[key].shape) == 1):
            x = d[key][0]
            y = d[key][1]
            #### Plot final cell position. MAY NOT TAKE RGB VALUES!!!
            plt.scatter(x, y, s = 5.0, c = d_attributes[key])

        else:
            print(key, " has no x,y points")

    #### Build plot frame, titles, and save data
    plt.ylim(-plot_y_extend/2, plot_y_extend/2)
    plt.xlim(-plot_x_extend/2, plot_x_extend/2)

    output_folder = ''
    snapshot = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
    snapshot = 'output' + f'{naming_index:08}'

    # Produce plot following the available options.
    if produce_for_panel == False:
        title_str = "History from image " + str(starting_index) + " to image " + str(endpoint) + "; " + title_str
        # %"Starting at frame {}, sample interval of {} for {} total samples".format(number_of_samples, sample_step_interval, number_of_samples)
        plt.title(title_str)
    else:
        fig.tight_layout()
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

    if output_plot is True:
        plt.savefig(output_folder + snapshot + '.png', dpi=256)
    if show_plot is True:
        plt.show()
    plt.close() # cell_tracker_movie.py:151: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).


def create_tracks_movie(data_path: str, save_path: str, save_name: str, start_file_index: int, end_file_index: int,
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
        j = i + 1 # this offsets the index so that we don't report that 0 samples have been taken, while stil producing an image.
        starting_index = j - max_number_of_samples

        #### Goes with "trail closing" block - not currently being used.
        projected_upper_sample_index = max_number_of_samples + starting_index
        max_samples_left = len(truncated_list_of_svgs) - j

        if i >= max_number_of_samples:
            plot_cell_tracks_for_movie(starting_index, 1, max_number_of_samples, output_plot, show_plot, i, produce_for_panel)
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


if __name__ == '__main__':
    # Execute only if run as a script

    if len(sys.argv) == 1:
        # Running with no arguments will make the script run every SVG with not stop to trail length

        create_tracks_movie('.', '', 'cell_tracks', 0, 10, 1, True, True)

    elif len(sys.argv) == 2:
        # Running with 1 argument sets the movie name and nothign else

        movie_name = sys.argv[1]
        create_tracks_movie('.', '', movie_name, 0, 10, 1, True, True)

    elif len(sys.argv) == 7:
        starting_file_index = int(sys.argv[1])
        end_file_index = int(sys.argv[2])
        cell_trail_length = int(sys.argv[3])  # length in time steps
        movie_name = sys.argv[4]
        INCLUDE_ALL_SVGs = bool(distutils.util.strtobool(sys.argv[5]))# bool(sys.argv[5])
        INCLUDE_FULL_HISTORY = bool(distutils.util.strtobool(sys.argv[6])) # bool(sys.argv[6])

        create_tracks_movie('.', '', movie_name, starting_file_index, end_file_index, cell_trail_length, INCLUDE_ALL_SVGs, INCLUDE_FULL_HISTORY)

    else:
        print('\nInput 0 arguments to process every available full and include full history and output movie with '
              'default name of cell_tracks.mp4')
        usage_str = "Usage: %s \n" % (sys.argv[0])
        print(usage_str)
        print('Input 1 argument (a string) to set movie name and process all files and full history')
        usage_str = "Usage: %s this_is_great_data\n" % (sys.argv[0])
        print(usage_str)

        print('Input 7 arguments to gain the most control')
        usage_str = "Usage: %s <start tracking index> <end file index> <history trail length> <movie name> " \
                    "<Bool: Include all SVGs> <Bool: Include full cell history> \n" % (sys.argv[0])
        print(usage_str)

        exit(1)
