# Plot the density of the ecm

from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sci
import math
from scipy import stats

sim_frame = 'output00001192' # The name of the simulation frame (eg. `initial`, `output00001192`).  This is not the file name
output_path = 'example_output' # The folder with all the output data

def load_ecm_data(file_name):
    mat = sci.loadmat(file_name)['ECM_Data']
    data = {"x":mat[0], "y":mat[1], "z":mat[2], "anisotropy":mat[3], "density":mat[4], "fiber_alignment_x":mat[5], "fiber_alignment_y":mat[6], "fiber_alignment_z":mat[7], "gradient_x":mat[8], "gradient_y":mat[9], "gradient_z":mat[10]}

    return data

# load data
mcds = pyMCDS(sim_frame + '.xml', output_path)
ecm_data = load_ecm_data(output_path + '/' + sim_frame + '_ECM.mat')

# Set our z plane and get our substrate values along it
z_val = 0.00
ecm_density = ecm_data['density']

# Get the 2D mesh for contour plotting
xy_len = int(math.sqrt(len(ecm_data['x'])))
xx, yy = np.array(ecm_data['x'][0:xy_len]), np.array(ecm_data['y'][0:xy_len])

# Make ecm density 3d
ecm_density = np.array(ecm_density)
ecm_density = ecm_density.reshape((xy_len, xy_len))

print(stats.describe(ecm_density))

# We want to be able to control the number of contour levels so we
num_levels = 21

# set up the figure area and add data layers
fig, ax = plt.subplots()
cs = ax.contourf(xx, yy, ecm_density, levels=num_levels)
ax.contour(xx, yy, ecm_density, color='black', levels = num_levels, linewidths=0.5)

# Now we need to add our color bar
cbar1 = fig.colorbar(cs, shrink=0.75)
cbar1.set_label('??')

# Let's put the time in to make these look nice
ax.set_aspect('equal')
ax.set_xlabel('x (micron)')
ax.set_ylabel('y (micron)')
ax.set_title('ECM density (??) at t = {:.1f} {:s}, z = {:.2f} {:s}'.format(mcds.get_time(), mcds.data['metadata']['time_units'], z_val, mcds.data['metadata']['spatial_units']))

plt.show()