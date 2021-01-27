import sys
import xml.etree.ElementTree as ET
import numpy as np
import glob
import matplotlib.pyplot as plt
import math, os, sys, re


def plot_cell_tracks(starting_index: int, sample_step_interval: int, number_of_samples: int, output_plot: bool,
                     show_plot: bool, naming_index: int):
    ####################################################################################################################
    ####################################            Usage example and input loading             ########################
    ####################################################################################################################

    # if (len(sys.argv) < 4):
    #     usage_str = "Usage: %s <start tracking index> <step interval for tracking> <# of samples to include>" % (sys.argv[0])
    #     print(usage_str)
    #     print("e.g.,")
    #     eg_str = "%s 0 1 10 indicates start at 0, go up by ones, and stop when you 10 samples" % (sys.argv[0])
    #     print(eg_str)
    #     exit(1)
    # else:
    #     starting_index = int(sys.argv[1])
    #     sample_step_interval = int(sys.argv[2])
    #     number_of_samples = int(sys.argv[3])
    #
    #
    # maxCount = starting_index

    d = {}  # dictionary to hold all (x,y) positions of cells

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

            if ('width' in child.attrib.keys()):
                #### Assumes square domains
                plot_spatial_length = float(child.attrib['width'])
                # print(plot_spatial_length)

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
                    d[child.attrib['id']] = np.vstack(
                        (d[child.attrib['id']], [float(circle.attrib['cx']), float(circle.attrib['cy'])]))
                else:
                    d[child.attrib['id']] = np.array([float(circle.attrib['cx']), float(circle.attrib['cy'])])
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

            #### Plot final cell position
            plt.scatter(x[-1], y[-1], s=3.5)

        elif (len(d[key].shape) == 1):
            x = d[key][0]
            y = d[key][1]

            plt.scatter(x, y, s=3.5)

        else:
            print(key, " has no x,y points")

    #### Build plot frame, titles, and save data

    ####### How do I get 1000 from the data itself???

    # xmlns: dc = "http://purl.org/dc/elements/1.1/"
    # xmlns: cc = "http://creativecommons.org/ns#"
    # xmlns: rdf = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    # xmlns: svg = "http://www.w3.org/2000/svg"
    # xmlns = "http://www.w3.org/2000/svg"
    # version = "1.1"
    # width = "1000"
    # height = "1070"
    plt.ylim(0, plot_spatial_length)
    plt.xlim(0, plot_spatial_length)

    title_str = "Starting at frame {}, sample interval of {} for {} total samples".format(starting_index,
                                                                                          sample_step_interval,
                                                                                          number_of_samples)
    plt.title(title_str)

    output_folder = ''
    snapshot = str(starting_index) + '_' + str(sample_step_interval) + '_' + str(number_of_samples)
    snapshot = 'output' + f'{naming_index:08}'

    #### Flags for output
    # output_plot = True
    # show_plot = True

    # Plot output
    if output_plot is True:
        plt.savefig(output_folder + snapshot + '.png')
    if show_plot is True:
        plt.show()
    # plt.close()


# if (len(sys.argv) < 2):
#     usage_str = "Usage: %s <start tracking index> <step interval for tracking> <# of samples to include>" % (sys.argv[0])
#     print(usage_str)
#     print("e.g.,")
#     eg_str = "%s 0 1 10 indicates start at 0, go up by ones, and stop when you 10 samples" % (sys.argv[0])
#     print(eg_str)
#     exit(1)
# else:
#     starting_index = int(sys.argv[1])
#     sample_step_interval = int(sys.argv[2])
#     number_of_samples = int(sys.argv[3])

# plot_cell_tracks(0, 1, 10, True, True)

def create_tracks_movie(data_path: str, save_path: str, save_name: str, start_file_index: int, end_file_index: int,
                        trail_length: int):
    """
    Generates the list of files in data_path, finds the ones with ECM data, makes plots from them, then outputs an
    ffmpeg generated movie to save_path, naming the movie save_name.

    This function requires ffmpeg be installed at the command line.

    :param data_path: Path to direcotry containing data
    :param save_path: Path to save generated image and movie to
    :param save_name: Save name for movie
    :return:
    """

    files = os.listdir(data_path)

    list_of_svgs = []

    for i in range(len(files)):
        if not re.search('snapshot(.*)\.svg', files[i]):
            continue

        # I feel like a dictionary could be used here, but I really need some ordering. A dict might be faster, but I don't
        # expect huge file lists. So I will just sort as I know how to do that ...

        list_of_svgs.append(files[i])

        # print(files[i])

    list_of_svgs.sort()

    truncated_list_of_svgs = []

    for i in range(len(list_of_svgs)):

        if i < start_file_index:
            continue

        if i >= end_file_index:
            continue

        truncated_list_of_svgs.append(list_of_svgs[i])

    # print(list_of_svgs)
    print(truncated_list_of_svgs)

    ALL_FILES = False

    if ALL_FILES:
        print('I went into the boolean')
        truncated_list_of_svgs = list_of_svgs

    max_number_of_samples = trail_length

    # Would be ideal to have a range - that would help with all the starting index, max samples left etc. Also, as written it isn't very flexible
    # would certainly be ideal to not call plot_cell_tracks every time, but instead store what is available. Could add a function that just
    # extracts the data from one SVG then appends it to exsisting data structure. could read all the desired data into Pandas DF
    # then write out images. Etc. Butas is, this is definitely reading the SVGs much to frequently.

    for i in range(len(truncated_list_of_svgs)):
        j = i + 1 # this offsets the index so that we don't report that 0 samples have been taken, while stil producing an image.
        print(j)
        starting_index = j - max_number_of_samples
        projected_upper_sample_index = max_number_of_samples + starting_index
        max_samples_left = len(truncated_list_of_svgs) - j
        print(len(truncated_list_of_svgs))

        if i >= max_number_of_samples:
            plot_cell_tracks(starting_index, 1, max_number_of_samples, True, True, i)
            print('middle')

        #### If one wanted to make the trails collapse into the last available location of the cell you would use something
        #### like this elif block
        # elif projected_upper_sample_index > len(list_of_svgs)-1:
        #     plot_cell_tracks(starting_index, 1, max_samples_left, True, True, i)
        #     print(max_samples_left)
        #     print('late')
        else:
            plot_cell_tracks(0, 1, j, True, True, i)
            print('early')

        # if i < max_number_of_samples:
        #     plot_cell_tracks(0, 1, j, True, True, i)
        #     print('early')
        #
        # #### If one wanted to make the trails collapse into the last available location of the cell you would use something
        # #### like this elif block
        # # elif projected_upper_sample_index > len(list_of_svgs)-1:
        # #     plot_cell_tracks(starting_index, 1, max_samples_left, True, True, i)
        # #     print(max_samples_left)
        # #     print('late')
        # else:
        #     plot_cell_tracks(starting_index, 1, max_number_of_samples, True, True, i)
        #     print('middle')
    # 11111
    # So close - need to get the file name right now ...
    number_frames = end_file_index - start_file_index

    string_of_interest = 'ffmpeg -start_number ' + str(
        start_file_index) + ' -y -framerate 6 -i ' + save_path + 'output%08d.png' + ' -frames:v ' + str(
        number_frames) + ' -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"'
    print(string_of_interest)
    os.system(
        'ffmpeg -start_number ' + str(
            start_file_index) + ' -y -framerate 6 -i ' + save_path + 'output%08d.png' + ' -frames:v ' + str(
            number_frames) + ' -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"')

    # https://superuser.com/questions/666860/clarification-for-ffmpeg-input-option-with-image-files-as-input
    # https://superuser.com/questions/734976/ffmpeg-limit-number-of-images-converted-to-video


if __name__ == '__main__':
    # Execute only if run as a script
    #

    print(len(sys.argv))
    if len(sys.argv) == 1:
        # Want this to be the equivalent of "run all SVGs present and make a movie out of all of them"
        # to do this - need a child of or an overload of "create_tracks_movie" that just runs everything

        create_tracks_movie('.', '', movie_name, starting_file_index, end_file_index, cell_trail_length)

    elif len(sys.argv) == 5:
        # usage_str = "Usage: %s <start tracking index> <stop tracking index> <# of samples to include>" % (sys.argv[0])
        # print(usage_str)
        # print("e.g.,")
        # eg_str = "%s 0 1 10 indicates start at 0, go up by ones, and stop when you 10 samples" % (sys.argv[0])
        # print(eg_str)
        starting_file_index = int(sys.argv[1])
        end_file_index = int(sys.argv[2])
        cell_trail_length = int(sys.argv[3])  # length in time steps
        movie_name = sys.argv[4]

        create_tracks_movie('.', '', movie_name, starting_file_index, end_file_index, cell_trail_length)
    else:
        exit(1)

    # create_tracks_movie('.', '', movie_name, starting_file_index, end_file_index, cell_trail_length)

# for fname in glob.glob('snapshot*.svg'):
#     print(fname)
#
#     start_index = int(fname[8:16])
#     print(start_index)
#
#     plot_cell_tracks(start_index, 1, 1, True, True)
# for fname in['snapshot00000000.svg', 'snapshot00000001.svg']:
# for fname in['snapshot00000000.svg']:
#   print(fname)
#   count += 1
#   if count > maxCount:
#     break
