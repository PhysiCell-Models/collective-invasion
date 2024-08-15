#!/bin/bash

# Script for making simple leader-follower movies and stills. Assumes data is already run!!!

# 0.1 um/min Output every 120 minutes
cp python_imaging/partial_history_2_level_contour_movie.py adh_10_replusion_25_speed_010/
cp python_imaging/partial_history_2_level_contour_still.py adh_10_replusion_25_speed_010/
cd adh_10_replusion_25_speed_010
python partial_history_2_level_contour_movie.py
cp two_contour_movie.mp4 ../1_images_for_paper/adh_10_replusion_25_speed_010.mp4
python partial_history_2_level_contour_still.py  24 24 # code requires two inputs - so just putting in the same number twice (only need one)
cp two_contour_still_24.png ../1_images_for_paper/discrete_remodeling_speed_010_discrete_model_figure_2880_minutes.png

cd ..

# 0.5 um/min Output every 30 minutes
cp python_imaging/partial_history_2_level_contour_movie.py adh_10_replusion_25_speed_050/
cp python_imaging/partial_history_2_level_contour_still.py adh_10_replusion_25_speed_050/
cd adh_10_replusion_25_speed_050
python partial_history_2_level_contour_movie.py
cp two_contour_movie.mp4 ../1_images_for_paper/adh_10_replusion_25_speed_050.mp4
python partial_history_2_level_contour_still.py  96 384
cp two_contour_still_96.png ../1_images_for_paper/discrete_remodeling_speed_050_discrete_model_figure_2880_minutes.png
cp two_contour_still_384.png ../1_images_for_paper/discrete_remodeling_speed_050_discrete_model_figure_11520_minutes.png

cd ..

# 0.8 um/min Output every 6 minutes
cp python_imaging/partial_history_2_level_contour_movie.py adh_10_replusion_25_speed_080/
cp python_imaging/partial_history_2_level_contour_still.py adh_10_replusion_25_speed_080/
cd adh_10_replusion_25_speed_080
python partial_history_2_level_contour_movie.py
cp two_contour_movie.mp4 ../1_images_for_paper/adh_10_replusion_25_speed_080.mp4
python partial_history_2_level_contour_still.py  150 150 # code requires two inputs - so just putting in the same number twice (only need one)
cp two_contour_still_150.png ../1_images_for_paper/discrete_remodeling_speed_080_discrete_model_figure_900_minutes.png



