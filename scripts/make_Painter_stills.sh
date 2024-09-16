#!/bin/bash

# Script for making simple test movies and images. Assumes data is already run!!!

cp python_imaging/simple_test_stills_cells_and_environment_Painter.py simple_tests_Painter_mixed/
cp python_imaging/simple_test_stills_cells_and_environment_Painter.py simple_tests_Painter_orthgonal/
cp python_imaging/simple_test_stills_cells_and_environment_Painter.py simple_tests_Painter_parallel/
cp python_imaging/simple_test_stills_cells_and_environment_Painter.py simple_tests_Painter_random/

# Simple Mixed
cd simple_tests_Painter_mixed
python simple_test_stills_cells_and_environment_Painter.py dont_make_inset
cp still_contour_and_cells.png ../1_images_for_paper/still_contour_and_cells_mix.png

# Simple test 1
cd ..

cd simple_tests_Painter_orthgonal
python simple_test_stills_cells_and_environment_Painter.py make_inset
cp still_contour_and_cells.png ../1_images_for_paper/still_contour_and_cells_orthogonal.png


# Simple test 2
cd ..

cd simple_tests_Painter_parallel
python simple_test_stills_cells_and_environment_Painter.py make_inset
cp still_contour_and_cells.png ../1_images_for_paper/still_contour_and_cells_parallel.png

# Simple test 3
cd ..

cd simple_tests_Painter_random
python simple_test_stills_cells_and_environment_Painter.py make_inset
cp still_contour_and_cells.png ../1_images_for_paper/still_contour_and_cells_random.png

cd ..

