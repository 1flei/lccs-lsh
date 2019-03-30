# lcsb
Longest Common Sub-sequence Bucketing algorithm


To build the project, use the following instructions:
mkdir build
cd build
cmake ..
make -j



To test on the Mnist dataset:
First compute the ground truth using
./lcsb -A ground_truth_angle -n 60000 -d 50 -q 1000 -D ../data/Mnist/Mnist.ds -Q ../data/Mnist/Mnist.q -O ../data/ts/Mnist.angle

Then run
