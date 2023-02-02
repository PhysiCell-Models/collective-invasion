import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *


mf = PhysiCellPlotter()

mf.produce_movie()