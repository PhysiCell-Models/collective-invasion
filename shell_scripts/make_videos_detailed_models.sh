#!/bin/bash

# call from PhysiCell root directory
# Leader follower

# makes videos straigh from the SVGs

cd leader_follower/
ffmpeg -r 24 -f image2 -i frame%04d.png -vcodec libx264 -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" -strict -2 -tune animation -crf 15 -acodec none leader_follower.mp4
cp leader_follower.mp4 ../1_images_for_paper/leader_follower_SVG.mp4

# Ductal invasion
cd ..

cd invasive_spheroid_output/
ffmpeg -r 24 -f image2 -i frame%04d.png -vcodec libx264 -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" -strict -2 -tune animation -crf 15 -acodec none invasion.mp4
cp invasion.mp4 ../1_images_for_paper/invasion_SVG.mp4

# Fibrosis
cd ..

cd fibrosis_test/
ffmpeg -r 24 -f image2 -i frame%04d.png -vcodec libx264 -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" -strict -2 -tune animation -crf 15 -acodec none fibrosis.mp4
cp fibrosis.mp4 ../1_images_for_paper/fibrosis_SVG.mp4



