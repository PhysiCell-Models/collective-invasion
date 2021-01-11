#
# cell_tracks.py - plot 2-D cell tracks associated with PhysiCell .svg files
#
# Usage:
#  python cell_tracks.py <max # of .svg frames>
# 
# Dependencies include matplotlib and numpy. We recommend installing the Anaconda Python3 distribution.
#
# Examples (run from directory containing the .svg files):
#  python cell_tracks.py 100
#
# Author: Randy Heiland
#
import sys
import xml.etree.ElementTree as ET
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
# from pyMCDS_ECM import *
# import matplotlib.gridspec as gridspec


#print(len(sys.argv))
if (len(sys.argv) < 4):
  usage_str = "Usage: %s <start tracking index> <step interval for tracking> <# frames/sample to include>" % (sys.argv[0])
  print(usage_str)
  print("e.g.,")
  eg_str = "%s 0 1 10" % (sys.argv[0])
  print(eg_str)
  exit(1)
else:
   starting_index = int(sys.argv[1])
   sample_step_interval = int(sys.argv[2])
   number_of_samples = int(sys.argv[3])


maxCount = starting_index

d={}   # dictionary to hold all (x,y) positions of cells

""" 
--- for example ---
In [141]: d['cell1599'][0:3]
Out[141]: 
array([[ 4900.  ,  4900.  ],
       [ 4934.17,  4487.91],
       [ 4960.75,  4148.02]])
"""

# mcds = pyMCDS('output00000000.xml', '.')
# xx, yy = mcds.get_2D_mesh()  # Mesh

#### Make a vector of file indices to get ####
endpoint = sample_step_interval*number_of_samples
file_indices = np.linspace(starting_index, endpoint, num=number_of_samples, endpoint=False)
print(file_indices)
count = 0
for fname in glob.glob('snapshot*.svg'):
# for fname in['snapshot00000000.svg', 'snapshot00000001.svg']:
# for fname in['snapshot00000000.svg']:
#   print(fname)
  count += 1
  if count > maxCount:
    break

#  print('\n---- ' + fname + ':')
  tree=ET.parse(fname)

#  print('--- root.tag, root.attrib ---')
  root=tree.getroot()
#  print('--- root.tag ---')
#  print(root.tag)
#  print('--- root.attrib ---')
#  print(root.attrib)


#  print('--- child.tag, child.attrib ---')
  numChildren = 0
  for child in root:
#  print(child.tag, child.attrib)
#    print('attrib=',child.attrib)
#  if (child.attrib['id'] == 'tissue'):
    if ('id' in child.attrib.keys()):
#      print('-------- found tissue!!')
      tissue_parent = child
      break

#  print('------ search tissue')
  for child in tissue_parent:
#    print('attrib=',child.attrib)
    if (child.attrib['id'] == 'cells'):
#      print('-------- found cells!!')
      cells_parent = child
      break
    numChildren += 1


  num_cells = 0
#  print('------ search cells')
  for child in cells_parent:
#    print(child.tag, child.attrib)
#    print('attrib=',child.attrib)
    for circle in child:
#      print('  --- cx,cy=',circle.attrib['cx'],circle.attrib['cy'])
      xval = float(circle.attrib['cx'])

      # should we test for bogus x,y locations??
      if (math.fabs(xval) > 10000.):
        print("xval=",xval)
        break
      yval = float(circle.attrib['cy'])
      if (math.fabs(yval) > 10000.):
        print("xval=",xval)
        break
      if (child.attrib['id'] in d.keys()):
        d[child.attrib['id']] = np.vstack((d[child.attrib['id']], [ float(circle.attrib['cx']), float(circle.attrib['cy']) ]))
      else:
        d[child.attrib['id']] = np.array( [ float(circle.attrib['cx']), float(circle.attrib['cy']) ])
      break
#    if (child.attrib['id'] == 'cells'):
#      print('-------- found cells!!')
#      tissue_child = child
    num_cells += 1
  print(fname,':  num_cells= ',num_cells)


fig = plt.figure(figsize=(8,8))
ax = fig.gca()
ax.set_aspect("equal")
#ax.set_xticks([])
#ax.set_yticks([]);
#ax.set_xlim(0, 8); ax.set_ylim(0, 8)

#print 'dir(fig)=',dir(fig)
#fig.set_figwidth(8)
#fig.set_figheight(8)

for key in d.keys():
  if (len(d[key].shape) == 2):
    x = d[key][:,0]
    y = d[key][:,1]
    plt.plot(x,y)
  else:
    print(key, " has no x,y points")
#    print(" d[",key,"].shape=", d[key].shape)
#    print(" d[",key,"].size=", d[key].size)
#    print( d[key])


title_str = "Starting at frame {}, sample intervale of {} for {} total samples".format(starting_index, sample_step_interval, number_of_samples)
plt.title(title_str)
plt.show()

# starting_index = int(sys.argv[1])
# sample_step_interval = int(sys.argv[2])
# number_of_samples = int(sys.argv[3])
