#!/usr/bin/env bash
./create_batch.py;
for input_file in input/param_file_{small,large}_*.in
do 
    echo $input_file
    ./run_batch.sh $input_file
done
./plot_data.py;
alert;

