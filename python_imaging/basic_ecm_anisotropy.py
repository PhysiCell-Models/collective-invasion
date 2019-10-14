# Plot the anisotropy of ecm

from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

sim_frame = 'output00001192' # The name of the simulation frame (eg. `initial`, `output00001192`).  This is not the file name
output_path = 'example_output' # The folder with all the output data

# load data
mcds = pyMCDS(sim_frame + '.xml', output_path)

# Set our z plane and get our substrate values along it
z_val = 0.00
plane_ecm = mcds.get_concentrations('ECM anisotropy', z_slice=z_val)

# Get the 2D mesh for contour plotting
xx, yy = mcds.get_mesh(flat=True)

# We want to be able to control the number of contour levels so we
num_levels = 21

# set up the figure area and add data layers
fig, ax = plt.subplots()
cs = ax.contourf(xx, yy, plane_ecm, levels=num_levels)
ax.contour(xx, yy, plane_ecm, color='black', levels = num_levels, linewidths=0.5)

# Now we need to add our color bar
cbar1 = fig.colorbar(cs, shrink=0.75)
cbar1.set_label('mmHg')

# Let's put the time in to make these look nice
ax.set_aspect('equal')
ax.set_xlabel('x (micron)')
ax.set_ylabel('y (micron)')
ax.set_title('oxygen (mmHg) at t = {:.1f} {:s}, z = {:.2f} {:s}'.format(mcds.get_time(), mcds.data['metadata']['time_units'], z_val, mcds.data['metadata']['spatial_units']))

plt.show()