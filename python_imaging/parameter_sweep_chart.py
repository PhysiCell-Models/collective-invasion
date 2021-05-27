import math, os, sys, re
import xml.etree.ElementTree as ET
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.colorbar as colorbar
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle
import distutils.util
import csv

# from pyMCDS_ECM import *
try:
    from pyMCDS_ECM import *
except ImportError:
    from pyMCDS import *


fig, ax = plt.subplots(figsize=(10, 8))

parameter_sweep_test_points = []

# just do it manually - you would have it done if you had done it manually ... Jesus.

# with open('parameter_sweep_points.csv', mode='r') as csv_file:
#     csv_reader = csv.DictReader(csv_file)
#     line_count = 0
#     for row in csv_reader:
#         if line_count == 0:
#             print(f'Column names are {", ".join(row)}')
#             line_count += 1
#         print(f'\t{row["Parameter Set number"]} works in the {row["Cell adhesion strength"]} department, and was born in {row["Cell speed"]}.')
#         line_count += 1
#     print(f'Processed {line_count} lines.')

with open('parameter_sweep_points.csv', mode='r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        print(type(row["Cell adhesion strength"]))
        print(type(float(row["Cell speed"])))
        ax.scatter(float(row["Cell adhesion strength"]), float(row["Cell speed"]), c='black')
        print(f'\t{row["Parameter Set number"]} works in the {row["Cell adhesion strength"]} department, and was born in {row["Cell speed"]}.')
        line_count += 1
    print(f'Processed {line_count} lines.')

ax.grid(alpha=0.5, which='major')
ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
plt.savefig('parameter_sweep_points.svg')

# try to make wider!!!

plt.show()