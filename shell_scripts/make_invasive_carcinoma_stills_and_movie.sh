#!/bin/bash

# Script for making invasive carcinoma movies and stills. Assumes data is already run!!!

# fibrosis 
cp python_imaging/partial_history_multilevel_contour_movie_invasive_carcinoma.py invasive_carcinoma_output/
cp python_imaging/partial_history_multilevel_contour_still_density.py invasive_carcinoma_output/
cd invasive_carcinoma_output
# show denisty instead of anisotropy
# python partial_history_multilevel_contour_movie_invasive_carcinoma.py
# cp multi_color_movie.mp4 ../1_images_for_paper/invasive_carcinoma.mp4

# doing two at a time because thats how the script is set up
python partial_history_multilevel_contour_still_density.py --snapshot_IDs 0 5 10 15 30 60 120 240 --contour_range 0.0 1.0
cp multi_contour_still_0.png ../1_images_for_paper/invasive_carcinoma_0.png
cp multi_contour_still_5.png ../1_images_for_paper/invasive_carcinoma_5.png
cp multi_contour_still_10.png ../1_images_for_paper/invasive_carcinoma_10.png
cp multi_contour_still_15.png ../1_images_for_paper/invasive_carcinoma_15.png
cp multi_contour_still_30.png ../1_images_for_paper/invasive_carcinoma_30.png
cp multi_contour_still_60.png ../1_images_for_paper/invasive_carcinoma_60.png
cp multi_contour_still_120.png ../1_images_for_paper/invasive_carcinoma_120.png
cp multi_contour_still_240.png ../1_images_for_paper/invasive_carcinoma_240.png

cp just_colorbar.png ../1_images_for_paper/invasive_carcinoma_colorbar.png