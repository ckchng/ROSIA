# ROSIA

The implementation of ROSIA: Rotation-Search-Based Star Identification Algorithm. ROSIA is built on top of the source code obtained from https://cs.adelaide.edu.au/~aparra/project/pcr/, which is an implementation of the paper "Parra Bustos, Alvaro, Tat-Jun Chin, and David Suter. "Fast rotation search with stereographic projections for 3d registration." Proceedings of the IEEE conference on computer vision and pattern recognition. 2014."

Prerequsites:  
CMAKE  
g++  
Python3  
numba (optional) < To speed up operations. Comment out line 3 and 23 in 'catalog_preprocessing.py' if this library is not installed.  

To install:    
mkdir build
cd build
cmake ..
make

To create specific onboard catalog:  
./create_onboard_dataset.sh

Hyperparamters:
1) mag_thres - include only stars that are observable within the sensor's sensitivity/
2) dist_thres - (in degrees) parameters to decide which stars are too close to be considered as separate stars. 
3) k - the number of closest distances to constraint each star, k = 2 is the triplet constraint as proposed in the paper.

Format of the pre-processed catalog:  
x y z Magnitude HIP index  CD1  CD2

Legends:  
x,y,z - the star vectors in the inertia's coordinates  
HIP index - the ID of stars in the HIPPARCOS catalog
CD1 - closest angular distance
CD2 - second closest angular distance

To create testing data:    
./create_onboard_dataset.sh

Format of the testing data ('scene_id.txt'):    
x y z Magnitude HIP index  

Legends:    
x,y,z - the star vectors in the camera's coordinates

Format of the testing data ('rot.txt')    
The 3x3 rotation matrix that rotates the testing star vectors into catalog's inertia's coordinates.  

Hyperparamters:
1) dist_noise_sigma - the sigma of the angular distance noise.
2) mag_noise_sigma - the sigma of the magnitude noise.
3) num_false_stars - the number of false stars in each testing instance.
4) num_scenes - number of scenes to be generated.
5) hip_catalog_filename - the directory of the original hipparcos catalog file.

To run ROSIA:  
./demo.sh 

Hyperparamters:  
1) dist_th - angular distance threshold (radians) to facilitate angular uncertainty of the star vectors.
2) mag_th - magnitude threshold (magnitude) to facilitate magnitude uncertainty of the stars.

Format of the output text files:  
lines 1 to N-3    
x y z HIP index    

lines N-3 to N    
The (estimated) 3x3 rotation matrix that rotates the testing star vectors into catalog's inertia's coordinates.  


To run the multi-pole Star-ID algorithm [1][2]. Note that the implementation doesn't utilise the k-vector technique, hence it is much slower than the actual algorithm. However, the ID results are not affected.   

./MPA.sh

Format of the output text files.    
x y z HIP index  

If you use this code, please cite us. Thanks!

@article{chng2023rosia,
  title={ROSIA: Rotation-Search-Based Star Identification Algorithm}, 
  author={Chng, Chee-Kheng and Bustos, {\'A}lvaro Parra and McCarthy, Benjamin and Chin, Tat-Jun}, 
  journal={IEEE Transactions on Aerospace and Electronic Systems}, 
  volume={59}, 
  number={5}, 
  pages={6469--6484}, 
  year={2023}, 
  publisher={IEEE} 
}

[1] Schiattarella, Vincenzo, Dario Spiller, and Fabio Curti. "A novel star identification technique robust to high presence of false objects: The Multi-Poles Algorithm." Advances in Space Research 59.8 (2017): 2133-2147.
[2] https://www.gaussteam.com/wordpress/wp-content/uploads/2018/02/IAA-AAS-CU-17-05-02-Schiattarella.pdf


