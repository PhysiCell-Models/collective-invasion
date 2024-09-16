#!/bin/bash

# Script for making simple test movies and images. Assumes data is already run!!!

# Simple test 0
cp python_imaging/simple_test_movies_cells_only.py simple_test0/
cp python_imaging/simple_test_stills_cells_and_environment_only.py simple_test0/
cd simple_test0
python simple_test_movies_cells_only.py
cp simple_test_cells.mp4 ../1_images_for_paper/simple_test0_straight_ECM.mp4
python simple_test_stills_cells_and_environment_only.py
cp still_initial_step.png ../1_images_for_paper/one_way_interactions_straight_ECM_no_cues_0.png
cp still_150.png ../1_images_for_paper/one_way_interactions_straight_ECM_no_cues_150.png
cp still_417.png ../1_images_for_paper/one_way_interactions_straight_ECM_no_cues_417.png

# Simple test 1
cd ..

cp python_imaging/simple_test_march_movie.py simple_test1/
cp python_imaging/simple_test_stills_march.py simple_test1/
cd simple_test1
# try making in studio: python simple_test_march_movie.py
# try making in studio: cp simple_test_cells.mp4 ../1_images_for_paper/simple_test1_march.mp4
python simple_test_stills_march.py
cp march_90.png ../1_images_for_paper/one_way_interactions_march_90.png
cp march_500.png ../1_images_for_paper/one_way_interactions_march_500.png
cp march_1200.png ../1_images_for_paper/one_way_interactions_march_1200.png


# Simple test 2
cd ..

cp python_imaging/simple_test_movies_cells_only.py simple_test2/
cp python_imaging/simple_test_stills_cells_and_environment_only.py simple_test2/
cd simple_test2
python simple_test_movies_cells_only.py
cp simple_test_cells.mp4 ../1_images_for_paper/simple_test2_random_1_D_circles.mp4 
python simple_test_stills_cells_and_environment_only.py 150 417
cp still_initial_step.png ../1_images_for_paper/one_way_interactions_circular_ECM_no_cues_0.png
cp still_150.png ../1_images_for_paper/one_way_interactions_circular_ECM_no_cues_150.png
cp still_417.png ../1_images_for_paper/one_way_interactions_circular_ECM_no_cues_417.png

# Simple test 3
cd ..

cp python_imaging/simple_test_movies_cells_only.py simple_test3/
cp python_imaging/simple_test_stills_cells_and_environment_only.py simple_test3/
cd simple_test3
python simple_test_movies_cells_only.py
cp simple_test_cells.mp4 ../1_images_for_paper/simple_test3_directed_circular_motion.mp4
python simple_test_stills_cells_and_environment_only.py 150 417
cp still_initial_step.png ../1_images_for_paper/one_way_interactions_circular_ECM_w_chemical_cue_0.png
cp still_150.png ../1_images_for_paper/one_way_interactions_circular_ECM_w_chemical_cue_150.png
cp still_417.png ../1_images_for_paper/one_way_interactions_circular_ECM_w_chemical_cue_417.png

# Simple test 4
cd ..

cp python_imaging/simple_test_movies_cells_only.py simple_test4/
cp python_imaging/simple_test_stills_cells_and_environment_only.py simple_test4/
cd simple_test4
python simple_test_movies_cells_only.py
cp simple_test_cells.mp4 ../1_images_for_paper/simple_test4_split_ECM.mp4
python simple_test_stills_cells_and_environment_only.py 60 150
cp still_initial_step.png ../1_images_for_paper/one_way_interactions_split_ECM_w_chemical_cue_0.png
cp still_60.png ../1_images_for_paper/one_way_interactions_split_ECM_w_chemical_cue_060.png
cp still_150.png ../1_images_for_paper/one_way_interactions_split_ECM_w_chemical_cue_150.png

