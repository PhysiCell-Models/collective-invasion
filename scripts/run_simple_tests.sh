#!/bin/bash
make -j4
./AMIGOS-invasion config/simple_test0_straight_ECM.xml
cp config/simple_test0_straight_ECM.xml simple_test0/ 
./AMIGOS-invasion config/simple_test1_cell_march.xml
cp config/simple_test1_cell_march.xml simple_test1/ 
./AMIGOS-invasion config/simple_test2_random_1_D_circles.xml
cp config/simple_test2_random_1_D_circles.xml simple_test2/
./AMIGOS-invasion config/simple_test3_directed_circular_motion.xml
cp config/simple_test3_directed_circular_motion.xml simple_test3/
./AMIGOS-invasion config/simple_test4_split_ECM.xml
cp config/simple_test4_split_ECM.xml simple_test4/

