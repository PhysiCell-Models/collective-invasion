#!/bin/bash
make -j4
# rm simple_tests_Painter_mixed/ simple_tests_Painter_orthogonal/ simple_tests_Painter_parallel/ simple_tests_Painter_random/
./AMIGOS-invasion config/simple_test_Painter_mixed.xml
# cp config/simple_test_Painter_mixed.xml simple_test0/ 
./AMIGOS-invasion config/simple_test_Painter_orthogonal.xml
# cp config/simple_test1_cell_march.xml simple_test1/ 
./AMIGOS-invasion config/simple_test_Painter_parallel.xml
# cp config/simple_test2_random_1_D_circles.xml simple_test2/
./AMIGOS-invasion config/simple_test_Painter_random.xml
# cp config/simple_test3_directed_circular_motion.xml simple_test3/

