#!/bin/sh                                                                                                                                                                                                                                                                


echo "== Compiling...  =="
make clean

make all

echo " "
echo "== Calculating parameter =="
./cal_param -n 60000 -m 7

if [ -d "index" ]; then
    rm -r index
fi
mkdir index

echo " "
echo "== Testing using different random seeds =="
for i in 1 2 3 4 5 6 7 8 9 10
do
    #echo " "
    echo "= Test $i ="
    ./srs -I -d 128 -i index/ -m 7 -n 10000 -s ../data/Mnist/Mnist.ds -y f -e $i
    
    #echo " "
    #echo "== Processing query workload... =="
    ./srs -Q -c 4 -g ../data/Mnist/Mnist.l2 -i index/ -k 1 -q ../data/Mnist/Mnist.q -t 6 -r 0.299203
done
