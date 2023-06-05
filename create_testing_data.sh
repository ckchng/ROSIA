#!/bin/bash

black_list_dir='./data/blacklist.txt'
out_dir='./data/'

dist_noise_sigma='2.4e-4'
mag_noise_sigma='0.3'
num_false_stars='2'

mag_threshold='6'
num_scenes='3'
hip_catalog_filename='./data/hip_main.dat'

python3 ./data/data_generator.py $black_list_dir $out_dir $dist_noise_sigma $mag_noise_sigma $num_false_stars $mag_threshold $num_scenes $hip_catalog_filename
