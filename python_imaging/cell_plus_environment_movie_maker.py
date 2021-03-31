from pyMCDS_ECM import *
import numpy as np

# Script REQUIRES ffmpeg to make movei!!!!!!!

######## If using on remote system, uncomment this line below to load correct matplotlib backend ################
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import math, os, sys, re
import scipy
import distutils.util


def print_stats(arr):
    """
    Produces relevant statistical output to screen given an array of any dimension. It flattens the in row-major style,
    the default np.flatten behavior.

    :param arr: any dimensional array, but it probably makes the most sense to be a 2-d array
    :return: Prints to termminal the array mean, quartiles, min, and max.
    """

    print("Mean: ", np.mean(arr.flatten()))
    print("Q2 quantile of arr : ", np.quantile(arr, .50))
    print("Q1 quantile of arr : ", np.quantile(arr, .25))
    print("Q3 quantile of arr : ", np.quantile(arr, .75))
    print("Min : ", arr.min())
    print("Max : ", arr.max())


def create_plot(snapshot, folder, output_folder='.', output_plot=True, show_plot=False):
    """
    Creates a plot as per instructions inside the function. As of 10.13.20 this was a plot of ECM-organoid simulations:
    a base layer of a contour plot of either the anisotropy or the oxygen, the cells in the smulation as a scatter plot,
    and finally the ECM orientation overlaid with a quiver plot.

    Parameters
    ----------
    snapshot :
        Base name of PhysiCell output files - eg 'output00000275' --> 'output' + '%08d'
    folder : str
        Path to input data
    output_folder : str
        Path for image output
    output_plot : bool
        True = image file will be made. Image output is required for movie production.
    show_plot : bool
        True = plot is displayed. Expected to be false for large batches.
    Returns
    -------
    Nothing :
        Produces a png image from the input PhysiCell data.
    """
    ####################################################################################################################
    ####################################            Load data                                   ########################
    ####################################################################################################################
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', folder)

    # loads and reads ECM data
    mcds.load_ecm(snapshot + '_ECM.mat', folder)

    # Get cell positions and attributes, microenvironment, and ECM data for plotting.

    # Cells
    cell_df = mcds.get_cell_df()

    #### Diffusion microenvironment
    xx, yy = mcds.get_2D_mesh()  # Mesh
    plane_oxy = mcds.get_concentrations('oxygen', 0.0)  # Oxyen (used for contour plot)

    #### ECM microenvironment
    xx_ecm, yy_ecm = mcds.get_2D_ECM_mesh()  # Mesh
    plane_anisotropy = mcds.get_ECM_field('anisotropy', 0.0)  # Anistropy (used for scaling and contour plot)
    # plane_anisotropy = micro # Used for contour plot

    ####################################################################################################################
    ####################################            Preprocessing                               ########################
    ####################################################################################################################

    #### Helper varialbes and functions ######

    # Number of contours (could include as a parameter)
    num_levels = 10  # 25 works well for ECM, 38 works well for oxygen

    # Make levels for contours
    levels_o2 = np.linspace(1e-14, 38, num_levels)
    # levels_ecm = np.linspace(1e-14, 1.0, num_levels)
    levels_ecm = np.linspace(0.90, 0.93, num_levels) # for the march environment - need to especially highlight small changes in anistoropy. 

    # Old function and scripting to scale and threshold anisotorpy values for later use in scaling lenght of ECM fibers
    # for visualization purposes.

    # micro = plane_anisotropy
    # micro_scaled = micro
    #
    # def curve(x):
    #     #return (V_max * x) / (K_M + x)
    #     return 0.5 if x > 0.5 else x

    # for i in range(len(micro)):
    #     for j in range(len(micro[i])):
    #         #micro_scaled[i][j] = 10 *  math.log10(micro[i][j] + 1) / math.log10(2)
    #         micro_scaled[i][j] = curve(micro[i][j])

    ##### Process data for plotting - weight fibers by anisotropy, mask out 0 anisotropy ECM units, get cell radii and types

    # Anisotropy strictly runs between 0 and 1. Element by element mulitplication produces weighted lengths between 0 - 1
    # for vizualization

    scaled_ECM_x = np.multiply(mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0], plane_anisotropy)
    scaled_ECM_y = np.multiply(mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0], plane_anisotropy)

    # if we want the arrows the same length instead
    ECM_x = mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0]
    ECM_y = mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0]

    # mask out zero vectors
    mask = plane_anisotropy > 0.0001

    # get unique cell types and radii
    cell_df['radius'] = (cell_df['total_volume'].values * 3 / (4 * np.pi)) ** (1 / 3)
    types = cell_df['cell_type'].unique()
    colors = ['yellow', 'blue']

    ####################################################################################################################
    ####################################            Plotting                                    ########################
    ####################################################################################################################

    # start plot and make correct size
    fig, ax = plt.subplots(figsize=(12, 12))
    plt.ylim(-500, 500)
    plt.xlim(-500, 500)

    # add contour layer
    # cs = plt.contourf(xx, yy, plane_oxy, cmap="Greens_r", levels=levels_o2)
    cs = plt.contourf(xx_ecm, yy_ecm, plane_anisotropy, cmap="YlGnBu", levels=levels_ecm)

    # Add cells layer
    for i, ct in enumerate(types):
        plot_df = cell_df[cell_df['cell_type'] == ct]
        for j in plot_df.index:
            circ = Circle((plot_df.loc[j, 'position_x'], plot_df.loc[j, 'position_y']),
                           radius=plot_df.loc[j, 'radius'], alpha=0.7, edgecolor='black')
            ax.add_artist(circ)

    # add quiver layer with scaled arrows ###
    # q = ax.quiver(xx_ecm[mask], yy_ecm[mask], scaled_ECM_x[mask], scaled_ECM_y[mask], pivot='middle', angles='xy', scale_units='inches', scale=2.0, headwidth=0,
    #               width=0.0015)  ## What is the deal with the line segment lengths shifting as the plots progress when I don't ue teh scaling??

    # add ECM orientation vectors unscaled by anistorpy ###
    plt.quiver(xx, yy, ECM_x, ECM_y,
    pivot='middle', angles='xy', scale_units='inches', scale=3.0, headwidth=0)

    # ax.axis('scaled') #used to be 'equal' https://stackoverflow.com/questions/45057647/difference-between-axisequal-and-axisscaled-in-matplotlib
    # This changes teh axis from -750,750 to ~-710,730. It looks better with scaled compared to axix, but either way it changes the plot limits

    # Labels and title (will need removed for journal - they will be added manually)
    ax.set_xlabel('x [micron]')
    ax.set_ylabel('y [micron]')
    fig.colorbar(cs, ax=ax)
    plt.title(snapshot)

    # Carefully place the command to make the plot square AFTER the color bar has been added.
    ax.axis('scaled')
    fig.tight_layout()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-500, 500)
    plt.xlim(-500, 500)

    # Plot output
    if output_plot is True:
        plt.savefig(output_folder + snapshot + '.png')
    if show_plot is True:
        plt.show()
    # plt.close()

def create_movie(data_path: str, save_path: str, save_name: str):
    """
    Generates the list of files in data_path, finds the ones with ECM data, makes plots from them, then outputs an
    ffmpeg generated movie to save_path, naming the movie save_name.

    This function requires ffmpeg be installed at the command line.

    :param data_path: Path to direcotry containing data
    :param save_path: Path to save generated image and movie to
    :param save_name: Save name for movie
    :return:
    """

    # generate list of files in the directory
    files = os.listdir(data_path)
    # files = list(filter(re.compile(r'output*ECM\.mat').search, files))

    # For all files in the directory, process only those with with 'ECM.mat' in the names. I am not sure why there is a
    # period at the beginning of the search pattern.
    for i in range(len(files)):
        if not re.search('.ECM\.mat', files[i]):
            continue

        # Sample call with meaningful variables:
        # create_plot('output00000275', output_folder='21_03_leader_follower_model_3_test/',output_plot=False, show_plot=False)
        create_plot(files[i].split('_')[0], data_path, output_folder=save_path, output_plot=True, show_plot=False)

    # make the movie - see ffmpeg documentation for more information

    # consider saving as jpegs - https://blender.stackexchange.com/questions/148231/what-image-format-encodes-the-fastest-or-at-least-faster-png-is-too-slow
    # consider compiling as movie instead of saving the files (all to increase processing speed) (then again, it was teh same speed)

    # consider not loading the unneeded data - and be sure to get rid of the unneeded fields!!!

    os.system(
        'ffmpeg -y -framerate 24 -i ' + save_path + 'output%08d.png -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"')


if __name__ == '__main__':
    # auto call the create movive function using the current directory as the data path and save path, and with teh given name.

    name_of_movie = sys.argv[1]

    create_movie('.', '', name_of_movie)