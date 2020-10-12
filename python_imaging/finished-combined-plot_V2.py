from pyMCDS_ECM import *
import numpy as np

######## If using on remote system, uncomment this line below to load correct matplotlib backend ################
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import math, os, sys, re
import scipy

def print_stats(arr):
    print("Mean: ", np.mean(arr.flatten()))
    print("Q2 quantile of arr : ", np.quantile(arr, .50)) 
    print("Q1 quantile of arr : ", np.quantile(arr, .25)) 
    print("Q3 quantile of arr : ", np.quantile(arr, .75))
    print("Min : ", arr.min())
    print("Max : ", arr.max())

def create_plot(snapshot, folder, output_folder = '.', output_plot = True, show_plot = False):
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', folder)
    # load ECM data
    mcds.load_ecm(snapshot + '_ECM.mat', folder)

    cell_df = mcds.get_cell_df()
    xx, yy = mcds.get_2D_mesh()
    xx_ecm, yy_ecm = mcds.get_2D_ECM_mesh()
    
    micro = mcds.get_ECM_field ('anisotropy', 0.0)

    # find levels for microenvironment
    plane_oxy = mcds.get_concentrations('oxygen', 0.0)
    plane_anisotropy = mcds.get_ECM_field('anisotropy', 0.0)
    num_levels = 25
    #levels = np.linspace(plane_oxy.min()+1e-14, plane_oxy.max(), num_levels)
    levels_o2 = np.linspace(1e-14, 38, num_levels)
    levels_ecm = np.linspace(1e-14, 1.0, num_levels)
    # arrow lengths depend on anisotropy
    micro_scaled = micro
    #print_stats(micro_scaled)
    #mean = np.mean(micro_scaled.flatten())
    V_max = 4
    #K_M = mean
    K_M = 0.4
    def curve(x):
        #return (V_max * x) / (K_M + x)
        return 0.5 if x > 0.5 else x

    # for i in range(len(micro)):
    #     for j in range(len(micro[i])):
    #         #micro_scaled[i][j] = 10 *  math.log10(micro[i][j] + 1) / math.log10(2)
    #         micro_scaled[i][j] = curve(micro[i][j])

    # micro_scaled = micro
    #print_stats(micro_scaled)

    # dy = mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0] * micro_scaled
    # dx = mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0] * micro_scaled

    dx = np.multiply(mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0], micro_scaled)
    dy = np.multiply(mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0], micro_scaled)

    #print(dx.shape)
    #print('dmag (min, max)', (np.sqrt(dx**2 + dy**2).min(), np.sqrt(dx**2 + dy**2).max()))

    # normalize lengths -- this needs some fiddling
    #dx = dx / dx.std()
    #dy = dy / dy.std()

    # if we want the arrows the same length instead
    dx_unscaled = mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0]
    dy_unscaled = mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0]

    # mask out zero vectors
    # mask = np.logical_or(dx > 1e-4, dy > 1e-4)
    mask = micro_scaled > 0.0001
    # get unique cell types and radii
    cell_df['radius'] = (cell_df['total_volume'].values * 3 / (4 * np.pi))**(1/3)
    types = cell_df['cell_type'].unique()
    colors = ['yellow', 'blue']

    fig, ax = plt.subplots(figsize=(12,12))
    y_limit_bottom_0, y_limit_top_0 = plt.ylim()
    x_limit_bottom_0, x_limit_top_0 = plt.xlim()
    # add contour layer
    # cs = plt.contourf(xx, yy, plane_oxy, cmap="Greens_r", levels=levels)
    cs = plt.contourf(xx_ecm, yy_ecm, plane_anisotropy, cmap="Reds", levels = levels_ecm)

    y_limit_bottom_1, y_limit_top_1 = plt.ylim(-800,800)
    x_limit_bottom_1, x_limit_top_1 = plt.xlim(-800,800)

    # Add cells layer
    for i, ct in enumerate(types):
        plot_df = cell_df[cell_df['cell_type'] == ct]
        for j in plot_df.index:
            circ = Circle((plot_df.loc[j, 'position_x'], plot_df.loc[j, 'position_y']),
                        color=colors[i], radius=plot_df.loc[j, 'radius'], alpha=0.7)
            ax.add_artist(circ)
    y_limit_bottom_2, y_limit_top_2 = plt.ylim()
    x_limit_bottom_2, x_limit_top_2 = plt.xlim()

    # add quiver layer with scaled arrows ###
    # plt.quiver(xx_ecm[mask], yy_ecm[mask], dx[mask], dy[mask], pivot='middle', angles='xy', units='width', headwidth=0, width=.0015)
    q = ax.quiver(xx_ecm[mask], yy_ecm[mask], dx[mask], dy[mask], pivot='middle', angles='uv', units='width', headwidth=0,
               width=.0015) ## What is the deal with the line segment lengths shifting as the plots progress when I don't ue teh scaling??
    # add unscaled arrows ###
    # plt.quiver(xx[mask], yy[mask], dx_unscaled[mask], dy_unscaled[mask], 
                # pivot='mid', angles='xy', headwidth=3)
    y_limit_bottom_3, y_limit_top_3 = plt.ylim()
    x_limit_bottom_3, x_limit_top_3 = plt.xlim()
    # ax.axis('scaled') #used to be 'equal' https://stackoverflow.com/questions/45057647/difference-between-axisequal-and-axisscaled-in-matplotlib
        # This changes teh axis from -750,750 to ~-710,730. It looks better with scaled compared to axix, but either way it changes the plot limits
    y_limit_bottom_4, y_limit_top_4 = plt.ylim()
    x_limit_bottom_4, x_limit_top_4 = plt.xlim()

    ax.set_xlabel('x [micron]')
    ax.set_ylabel('y [micron]')
    fig.colorbar(cs, ax=ax) # --> this makes the plot super whonky and squishes our physicially realistic plot. Have to plot differnently.
    plt.title(snapshot)
    ax.axis('scaled')
    y_limit_bottom_4, y_limit_top_4 = plt.ylim(-800,800)
    x_limit_bottom_4, x_limit_top_4 = plt.xlim(-800,800)
    if output_plot is True:
        plt.savefig(output_folder + snapshot + '.png')
    if show_plot is True:
        plt.show()
    plt.close()

def create_gif(data_path, save_path, save_name):
    files = os.listdir(data_path)
    #files = list(filter(re.compile(r'output*ECM\.mat').search, files))
    for i in range(len(files)):
        if not re.search('.ECM\.mat', files[i]):
            continue

        create_plot(files[i].split('_')[0], data_path, output_folder=save_path, output_plot=True, show_plot=False)
        #exit()
    os.system('ffmpeg -y -framerate 24 -i ' + save_path + 'output%08d.png -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"')

if __name__ == '__main__':
    #create_plot('output00000275', output_folder='21_03_leader_follower_model_3_test/',output_plot=False, show_plot=True)
    create_gif('.', '', 'anisotropy')
    