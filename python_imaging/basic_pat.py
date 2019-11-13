from pyMCDS import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# load cell and microenvironment data
mcds = pyMCDS('output00000275.xml', 'Archive')
# load ECM data
mcds.load_ecm('output00000275_ECM.mat', 'Archive')

cell_df = mcds.get_cell_df()
xx, yy = mcds.get_2D_mesh()
micro = mcds.get_concentrations('ECM anisotropy', 0.0)

# find levels for microenvironment
num_levels = 20
levels = np.linspace(micro.min()+1e-14, micro.max(), num_levels)

# arrow lengths depend on anisotropy
dy = mcds.data['ecm']['y_vec'][:, :, 0] * micro
dx = mcds.data['ecm']['x_vec'][:, :, 0] * micro
print(dx.shape)
print('dmag (min, max)', (np.sqrt(dx**2 + dy**2).min(), np.sqrt(dx**2 + dy**2).max()))

# normalize lengths -- this needs some fiddling
dx = dx / dx.std()
dy = dy / dy.std()

# if we want the arrows the same length instead
dx_unscaled = mcds.data['ecm']['x_vec'][:, :, 0]
dy_unscaled = mcds.data['ecm']['y_vec'][:, :, 0]

# mask out zero vectors
mask = np.logical_or(dx > 1e-4, dy > 1e-4)

# get unique cell types and radii
cell_df['radius'] = (cell_df['total_volume'].values * 3 / (4 * np.pi))**(1/3)
types = cell_df['cell_type'].unique()
colors = ['yellow', 'blue']

fig, ax = plt.subplots(figsize=(12, 10))

# add contour layer
cs = plt.contourf(xx, yy, micro, cmap='Reds', levels=levels)

# Add cells layer
for i, ct in enumerate(types):
    plot_df = cell_df[cell_df['cell_type'] == ct]
    for j in plot_df.index:
        circ = Circle((plot_df.loc[j, 'position_x'], plot_df.loc[j, 'position_y']),
                       color=colors[i], radius=plot_df.loc[j, 'radius'], alpha=0.7)
        ax.add_artist(circ)

# add quiver layer with scaled arrows ###
plt.quiver(xx[mask], yy[mask], dx[mask], dy[mask], pivot='mid', angles='xy', scale=200,
                units='width', headwidth=2)

# add unscaled arrows ###
# plt.quiver(xx[mask], yy[mask], dx_unscaled[mask], dy_unscaled[mask], 
            # pivot='mid', angles='xy', headwidth=3)

ax.axis('equal')
ax.set_xlabel('x [micron]')
ax.set_ylabel('y [micron]')
fig.colorbar(cs, ax=ax)

plt.show()
#plt.savefig('vector.png')