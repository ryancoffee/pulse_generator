#!/bin/bash
source set_vars;
for i in 0 2 3 .. 10
do
./generate_dataset;
sleep 40;
done

