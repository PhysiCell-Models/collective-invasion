#!/bin/bash
make -j4
./AMIGOS-invasion config/leader_follower_instant_speed_010.xml
# cp config/leader_follower_instant_speed_010.xml adh_10_replusion_25_speed_010/ 

./AMIGOS-invasion config/leader_follower_instant_speed_050.xml
# cp config/leader_follower_instant_speed_050.xml adh_10_replusion_25_speed_050/ 

./AMIGOS-invasion config/leader_follower_instant_speed_080.xml
# cp config/leader_follower_instant_speed_080.xml adh_10_replusion_25_speed_080/