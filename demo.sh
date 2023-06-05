#!/bin/bash


dist_th='4.8e-4'
mag_th='0.6'

dataset_dir="./data/prepro_catalog.txt"

in_dir="./data/"
out_dir="./output/"

i='0'

taskset -c 1 ./build/ROSIA_c $dist_th $mag_th $dataset_dir $out_dir $in_dir $i



