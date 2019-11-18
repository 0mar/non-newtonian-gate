./create_batch.py;
for input in input/param_file_{small,medium,large}_*.in
do ./run_batch.sh $input;
done;
./plot_data.py;
alert;

