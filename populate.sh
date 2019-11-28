#!/usr/bin/env bash
# This script is called by make with two arguments: source directory and build directory (in that order)
# Because it is a pain to spell out the copy command in cmake for each file.
if [[ "$#" -ne 2 ]]
then
    echo "Illegal number of parameters"
    exit 1
fi
SOURCE="$1"
BUILD="$2"
cp "$SOURCE"/create_*.py "$SOURCE"/params_*.json "$SOURCE"/run_batch.sh "$SOURCE"/plot_data.py "$SOURCE"/visualise.py "$SOURCE"/*_full.sh "$BUILD"
mkdir -p single_channel_data/plots double_channel_data/plots
