#!/bin/bash
make 
echo "N,calculated_pi,error,time" > serial.dat

for i in $(seq 1 1 2000); do
    echo "Running with input $i"
    ./calcpi "$i" >> serial.dat
done
make clean
python3 serial.py
