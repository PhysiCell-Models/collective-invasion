#!/bin/bash

# Script for making fibrosis movies and stills. Assumes data is already run!!!

# fibrosis 
cp python_imaging/partial_history_multilevel_contour_movie_fibrosis.py fibrosis_test/
cp python_imaging/partial_history_multilevel_contour_movie_fibrosis.py fibrosis_test/
cd fibrosis_test
# show denisty instead of anisotropy
# python partial_history_multilevel_contour_movie_fibrosis.py
# cp multi_color_movie.mp4 ../1_images_for_paper/fibrosis.mp4

# doing two at a time because thats how the script is set up
python partial_history_multilevel_contour_still_density.py --snapshot_IDs 0 30 60 120 240 360 --contour_range 0.5 1.0
cp multi_contour_still_0.png ../1_images_for_paper/fibrosis_0.png
cp multi_contour_still_30.png ../1_images_for_paper/fibrosis_30.png

# python partial_history_multilevel_contour_still_density.py 60 120 0.5 1.0
cp multi_contour_still_60.png ../1_images_for_paper/fibrosis_60.png
cp multi_contour_still_120.png ../1_images_for_paper/fibrosis_120.png

# python partial_history_multilevel_contour_still_density.py 240 360 0.5 1.0
cp multi_contour_still_240.png ../1_images_for_paper/fibrosis_240.png
cp multi_contour_still_360.png ../1_images_for_paper/fibrosis_360.png

cp just_colorbar.png ../1_images_for_paper/fibrosis_colorbar.png
