#!/bin/bash
make -j4
./AMIGOS-invasion config/writing_only.xml
cp config/writing_only.xml adh_0_repulsion_0_speed_10_no_reading/ 

./AMIGOS-invasion config/reading_only.xml
cp config/reading_only.xml adh_0_repulsion_0_speed_10_no_writing/ 

./AMIGOS-invasion config/writing_and_reading.xml
cp config/writing_and_reading.xml adh_0_repulsion_0_speed_10/

