# ROSIA

The implementation of ROSIA: Rotation-Search-Based Star Identification Algorithm

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

Format of the pre-processed catalog.
1 2 3     4         5       6    7
x y z magnitude HIP index  CD1  CD2

x,y,z - the star vectors in the inertia's coordinates
CD1 - closest angular distance
CD2 - second closest angular distance

To create testing data:
./create_onboard_dataset.sh

Format of the testing data ('scene_id.txt')
1 2 3     4        5     
x y z magnitude HIP index

x,y,z - the star vectors in the camera's coordinates

('rot.txt')
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

Format of the output text files.

lines 1 to N-3
1 2 3     4 
x y z HIP index

lines N-3 to N
The 3x3 rotation matrix that rotates the testing star vectors into catalog's inertia's coordinates.



