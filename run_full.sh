./create_batch.py;
for input in param_file_{small,medium,large}_*.in
do ./run_batch.sh $input;
done;
./plot_data.py;
alert;

