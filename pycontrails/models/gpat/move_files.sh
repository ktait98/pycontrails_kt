#!/bin/bash

# Define source and destination directories
SOURCE_DIR="Kieran2024/outputs/"
DEST_DIR1="Kieran2024/outputs/id_vs_plume/"
DEST_DIR2="Kieran2024/outputs/mc_v/"
DEST_DIR3="Kieran2024/outputs/sensitivity/"
DEST_DIR4="Kieran2024/outputs/sensitivity/sensitivity_locs/"

dir=$1

# Loop through directories starting with 'sensitivity' in the source directory
if [ $dir == "id_vs_plume" ]; then
  for dir in "$SOURCE_DIR"id_vs_plume_*; 
    do
      # Rsync the directory to the destination and remove source files
      rsync -av --remove-source-files --progress "$dir" "$DEST_DIR1"
    done
  echo "$dir done"
elif [ $dir == "mc_v" ]; then
  for dir in "$SOURCE_DIR"/mc_v*; 
    do
      rsync -av --remove-source-files --progress "$dir" "$DEST_DIR2"
    done
  echo "$dir done"
elif [ $dir == "sensitivity" ]; then
  for dir in "$SOURCE_DIR"/sensitivity*; 
    do
      rsync -av --remove-source-files --progress "$dir" "$DEST_DIR3"
    done
  echo "$dir done"
elif [ $dir == "sensitivity_locs" ]; then
  for dir in "$SOURCE_DIR"/sensitivity/sensitivity_locs*; 
    do
      rsync -av --remove-source-files --progress "$dir" "$DEST_DIR4"
    done
  echo "$dir done"
fi