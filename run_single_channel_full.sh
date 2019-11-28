#!/usr/bin/env bash
./create_single_channel_batch.py;
for input_file in single_channel_data/param_file_{small,large}_*.in
do 
    echo $input_file
    ./single_channel_batch.sh $input_file
done
./plot_data.py;
alert;

