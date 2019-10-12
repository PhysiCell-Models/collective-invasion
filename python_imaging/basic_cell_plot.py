# Plot the the cells given a file

from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

file_name = 'output00001192.xml'
output_path = 'example_output'

# load data
mcds = pyMCDS(file_name, output_path)

# Set our z plane and get our substrate values along it
z_val = 0.00
plane_oxy = mcds.get_concentrations('oxygen', z_slice=z_val)

# Get the 2D mesh for contour plotting
xx, yy = mcds.get_mesh(flat=True)
 
# We want to be able to control the number of contour levels so we
# need to do a little set up
num_levels = 21
 
# get our cells data and figure out which cells are in the plane
cell_df = mcds.get_cell_df()
ds = mcds.get_mesh_spacing()
inside_plane = (cell_df['position_z'] < z_val + ds) & (cell_df['position_z'] > z_val - ds)
plane_cells = cell_df[inside_plane]
 
# We're going to plot two types of cells and we want it to look nice
colors = ['blue', 'yellow']
sizes = [20, 8]
labels = ['follower', 'leader']

# set up the figure area and add microenvironment layer
fig, ax = plt.subplots()

# get our cells of interest
followers = plane_cells[plane_cells['cell_type'] == 2] # Followers are type 2
leaders = plane_cells[plane_cells['cell_type'] == 1] # Leaders are type 1

# plot the cell layer
for i, plot_cells in enumerate((followers, leaders)):
    ax.scatter(plot_cells['position_x'].values, 
            plot_cells['position_y'].values, 
            facecolor='none', 
            edgecolors=colors[i],
            alpha=0.6,
            s=sizes[i],
            label=labels[i])
 
 
# Let's put the time in to make these look nice
ax.set_aspect('equal')
ax.set_xlabel('x (micron)')
ax.set_ylabel('y (micron)')
ax.set_title('z = {:.2f} {:s}'.format(z_val, mcds.data['metadata']['spatial_units']))
ax.legend(loc='upper right')
 
plt.show()
