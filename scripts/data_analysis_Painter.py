import numpy as np
import csv


try:
    from pyMCDS_ECM import *
except ImportError:
    from pyMCDS import *


## Adapting the code from the image processing code - create_density_histogram function. 
## I don't need the histogram - just the 95th percentile value - converted back from bin number to position.

def create_density_histogram(mcds):
    
    # inherited from the image processing code
    num_bins = 40

    y_positions = mcds.data['discrete_cells']['position_y']

    # Create a 1D histogram of the y-positions, ensuring the range spans -500 to 500
    hist, bin_edges = np.histogram(y_positions, bins=num_bins, range=(-500, 500)) ## not automatic

    # Print the histogram counts
    print("Histogram counts for each bin:")
    for i, count in enumerate(hist):
        print(f"Bin {i+1}: {count} counts")

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Calculate the cumulative sum and find the 95th percentile value
    cumulative_counts = np.cumsum(hist)
    
    total_count = cumulative_counts[-1]
    percentile_95 = np.searchsorted(cumulative_counts, 0.95 * total_count)

    print(f"95th percentile bin number: {percentile_95}")
    print(f"95th percentile center of bin in um: {bin_centers[percentile_95]}")
    # to doulbe check - print bin centers then bin edges
    # print(bin_centers)
    # print(bin_edges)

    return bin_centers[percentile_95]

percentile_centers = []
for n in range(1, 21):
    # folder name needs padded with 1 zero
    if n < 10:
        folder_name = "0" + str(n)
    else:
        folder_name = str(n)

    #     Parameters
    # ----------
    # xml_name: str
    #     String containing the name of the xml file without the path
    # output_path: str, optional
    #     String containing the path (relative or absolute) to the directory
    #     where PhysiCell output files are stored (default= ".")

    mcds = pyMCDS('output00000005.xml', folder_name)


    percentile_centers.append(create_density_histogram(mcds))

# add in value for the extra run (the run used for detail in the paper)

percentile_centers.append(0.0) # placeholder for the extra run - these were determined by outputing the value from the image processing code. so just make the image 
                            # processing code output the value and then copy it here.   

print(percentile_centers)
print('Mean of 95th percentile centers:')
print(np.mean(percentile_centers))

print('Standard deviation of 95th percentile centers:')
print(np.std(percentile_centers))

np.savetxt("95th_percentile_center.csv", percentile_centers, delimiter=",", fmt="%s") # rename as needed