#!/bin/bash
make 
echo "N,calculated_pi,error,time" > serial.dat
for i in 10 100 1000 10000 100000 1000000 10000000; do
    echo "Running with input $i"
    ./calcpi "$i" >> serial.dat
done
make clean
