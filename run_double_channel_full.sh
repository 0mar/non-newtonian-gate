#!/usr/bin/env bash
./create_double_channel_batch.py;
for input_file in double_channel_data/*.in
do 
    ./double_channel_batch.sh $input_file
done
#./plot_data.py;
alert;

