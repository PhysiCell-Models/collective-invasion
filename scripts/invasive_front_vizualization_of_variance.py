import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV files into DataFrames
df_parallel = pd.read_csv('../Painter_data/95th_percentile_center_parallel.csv', header=None, names=['Position'])
df_perpendicular = pd.read_csv('../Painter_data/95th_percentile_center_orthogonal.csv', header=None, names=['Position'])
df_random = pd.read_csv('../Painter_data/95th_percentile_center_random.csv', header=None, names=['Position'])

# Add a 'Category' column to each DataFrame to differentiate the data
df_parallel['Category'] = 'Parallel'
df_perpendicular['Category'] = 'Perpendicular'
df_random['Category'] = 'Random'

# Combine the DataFrames into a single DataFrame
df_combined = pd.concat([df_random, df_parallel, df_perpendicular])

# Assuming the positional data is in a column called 'Position' in the CSV files
# If the column name is different, replace 'Position' with the actual name
sns.set(style="whitegrid")

# Create the beeswarm plot to show variations
fig, ax = plt.subplots(figsize=(5, 6))
# ax = sns.stripplot(x='Category', y='Position', data=df_combined, size=8)

# # plot the mean line
# sns.boxplot(showmeans=True,
#             meanline=True,
#             meanprops={'color': 'k'},
#             medianprops={'visible': False},
#             whiskerprops={'visible': False},
#             zorder=10,
#             x="Category",
#             y="Position",
#             data=df_combined,
#             showfliers=False,
#             showbox=False,
#             showcaps=False,
#             ax=ax)



#  Boxplot to show distribution statistics, mean is shown as a green triangle
sns.boxplot(x='Category', y='Position', data=df_combined, showmeans=True, meanline=True, width=0.3, 
            # meanprops={"marker":".", "markerfacecolor":"red", "markeredgecolor":"red", "markersize":"5"},
            boxprops=dict(alpha=0.7), ax=ax)

sns.pointplot(x='Category', y='Position', data=df_combined, estimator=np.mean, join=False, color='red', markers='o', scale=0.75, ax=ax)

# Jitter plot (strip plot) overlayed on the boxplot
sns.stripplot(x='Category', y='Position', data=df_combined, jitter=True, size=5, color='black', alpha=0.5, ax=ax)

plt.tick_params(axis='x', which='major', labelsize=14)
plt.tick_params(axis='y', which='major', labelsize=14)

# Customize the plot

# Set y-axis intervals to 125
plt.yticks(range(-500, 501, 125))

# Add grid lines
ax.set_axisbelow(True) # https://stackoverflow.com/questions/1726391/matplotlib-draw-grid-lines-behind-other-graph-elements # Not workign - but overall fig gets job done. 

# plt.rc('axes', axisbelow=True)
ax.yaxis.grid(linestyle='--', linewidth=0.7)

plt.ylabel('Position (µm)', fontsize=16)
plt.xlabel('', fontsize=2)
plt.ylim(-500, 500)  # Set limits based on your preference
plt.tight_layout()
# Save figure
plt.savefig('painter_extent_of_invasive_front.png', dpi=300)

# Show the plot
plt.show()
