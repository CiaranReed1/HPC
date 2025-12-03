#!/bin/bash
make 
echo "N,calculated_pi,error,time" > serial.dat
for i in $(seq 10 1 5000); do
    echo "Running with input $i"
    ./calcpi "$i" >> serial.dat
done
make clean
python3 serial.py