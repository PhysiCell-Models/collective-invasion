#!/bin/bash

# Script for making simple leader-follower movies and stills. Assumes data is already run!!!

# Writing only (no reading)
cp python_imaging/partial_history_2_level_contour_movie.py adh_0_repulsion_0_speed_10_no_reading/
cp python_imaging/partial_history_2_level_contour_still.py adh_0_repulsion_0_speed_10_no_reading/
cd adh_0_repulsion_0_speed_10_no_reading
python partial_history_2_level_contour_movie.py
cp two_contour_movie.mp4 ../1_images_for_paper/adh_0_repulsion_0_speed_10_no_reading.mp4
python partial_history_2_level_contour_still.py 60 160
cp two_contour_still_60.png ../1_images_for_paper/adh_0_repulsion_0_speed_10_no_reading_60.png
cp two_contour_still_160.png ../1_images_for_paper/adh_0_repulsion_0_speed_10_no_reading_160.png

cd ..

# Reading only (no writing)
cp python_imaging/partial_history_2_level_contour_movie.py adh_0_repulsion_0_speed_10_no_writing/
cp python_imaging/partial_history_2_level_contour_still.py adh_0_repulsion_0_speed_10_no_writing/
cd adh_0_repulsion_0_speed_10_no_writing
python partial_history_2_level_contour_movie.py
cp two_contour_movie.mp4 ../1_images_for_paper/adh_0_repulsion_0_speed_10_no_writing.mp4
python partial_history_2_level_contour_still.py 60 160
cp two_contour_still_60.png ../1_images_for_paper/adh_0_repulsion_0_speed_10_no_writing_60.png
cp two_contour_still_160.png ../1_images_for_paper/adh_0_repulsion_0_speed_10_no_writing_160.png

cd ..

# Both reading and writing
cp python_imaging/partial_history_2_level_contour_movie.py adh_0_repulsion_0_speed_10/
cp python_imaging/partial_history_2_level_contour_still.py adh_0_repulsion_0_speed_10/
cd adh_0_repulsion_0_speed_10
python partial_history_2_level_contour_movie.py
cp two_contour_movie.mp4 ../1_images_for_paper/adh_0_repulsion_0_speed_10.mp4
python partial_history_2_level_contour_still.py 60 160
cp two_contour_still_60.png ../1_images_for_paper/adh_0_repulsion_0_speed_10_60.png
cp two_contour_still_160.png ../1_images_for_paper/adh_0_repulsion_0_speed_10_160.png



