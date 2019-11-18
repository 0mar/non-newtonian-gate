#!/usr/bin/env bash
# This script is called by make with two arguments: source directory and build directory (in that order)
# Because it is a pain to spell out the copy command in cmake for each file.
SOURCE="$1"
BUILD="$2"
echo $SOURCE
echo $BUILD
cp "$SOURCE"/create_batch.py "$SOURCE"/params.json "$SOURCE"/run_batch.sh "$SOURCE"/plot_data.py "$SOURCE"/visualise.py "$SOURCE"/run_full.sh "$BUILD"
mkdir -p input/plots
