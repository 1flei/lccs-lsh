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

 ./lcsb -A srp_lcs -n 60000 -q 100 -d 50 -K 512 -D ../data/Mnist/Mnist.ds -Q ../data/Mnist/Mnist.q -G ../data/ts/Mnist.angle -O ./srp_lcslsh.out
