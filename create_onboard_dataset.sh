#!/bin/bash

catalog_path='./data/HIP_full_cat.txt'
blacklist_path='./data/blacklist.txt'
prepro_cat_path='./data/prepro_catalog.txt'
mag_thres='6'
dist_thres='0.2' # (degree) to remove stars that are too close. 0.2 is equal to 1 pixel to the default setting (see code)
k='2' # the number of closest distances to constraint each star, k = 2 is the triplet constraint as proposed in the paper.

python3 ./data/catalog_preprocessing.py $catalog_path $blacklist_path $prepro_cat_path $mag_thres $dist_thres $k
