#!/bin/bash

# Script for making leader-follower movies and stills. Assumes data is already run!!!

# leader follower
cp python_imaging/partial_history_multilevel_contour_movie.py leader_follower/
cp python_imaging/partial_history_multilevel_contour_still.py leader_follower/
cd leader_follower
python partial_history_multilevel_contour_movie.py
cp multi_color_movie.mp4 ../1_images_for_paper/leader_follower.mp4
python partial_history_multilevel_contour_still.py 480 1920
cp multi_contour_still_480.png ../1_images_for_paper/continous_remodeling_speed_050_2880_minutes.png
cp multi_contour_still_1920.png ../1_images_for_paper/continous_remodeling_speed_050_11520_minutes.png

cd ..
cp python_imaging/partial_history_multilevel_contour_movie.py leader_follower_decreased_remodeling/
cp python_imaging/partial_history_multilevel_contour_still.py leader_follower_decreased_remodeling/
cd leader_follower_decreased_remodeling
python partial_history_multilevel_contour_movie.py
cp multi_color_movie.mp4 ../1_images_for_paper/leader_follower_decreased_remodeling.mp4
python partial_history_multilevel_contour_still.py 480 1920
cp multi_contour_still_480.png ../1_images_for_paper/continous_remodeling_decreased_remodeling_2880_minutes.png
cp multi_contour_still_1920.png ../1_images_for_paper/continous_remodeling_decreased_remodeling_11520_minutes.png








