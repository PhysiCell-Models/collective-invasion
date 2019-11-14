# Python Imaging

In an effort to create some better imaging to view cells, ECM components and other things simultaneously, these scripts have been written.


## Scripts

* `basic_cell_plot.py`: create a plot with Cells plotted.
* `basic_oxygen_density.py`: Plot the density of oxygen.  
* `basic_ecm_density.py`: Plot the ECM density.
* `basic_ecm_anisotropy.py`: Plot the ECM anisotropy.

## ToDo
* Read in data
* ECM anisotropy plot
* ECM alignment plot
* Overlap all elements / modularity
* Panels

## Other stuff

### Helpful function to load in ECM matlab data
```
    import scipy.io as sci

    sim_frame = 'output00001192' # The name of the simulation frame (eg. `initial`, `output00001192`).  This is not the file name
    output_path = 'example_output' # The folder with all the output data

    def load_ecm_data(file_name):
        mat = sci.loadmat(file_name)['ECM_Data']
        data = {"x":mat[0], "y":mat[1], "z":mat[2], "anisotropy":mat[3], "density":mat[4], "fiber_alignment_x":mat[5], "fiber_alignment_y":mat[6], "fiber_alignment_z":mat[7], "gradient_x":mat[8], "gradient_y":mat[9], "gradient_z":mat[10]}

        return data

    # load data
    ecm_data = load_ecm_data(output_path + '/' + sim_frame + '_ECM.mat')
```