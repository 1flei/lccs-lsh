#!/bin/bash
make
rm *.o

# ------------------------------------------------------------------------------
#  Public Parameters
# ------------------------------------------------------------------------------
dname=Gist
n=1000000
d=960
qn=1000
dPath=../data/${dname}/${dname}
oFolder=../results/${dname}/

# ------------------------------------------------------------------------------
#  Ground-Truth
# ------------------------------------------------------------------------------
# ./lcsb -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
#   -ts ${dPath}.mcs

# ------------------------------------------------------------------------------
#  Algorithms for Maximum Cosine Similarity Search (MCSS)
#  01 - linear scan
#  02 - m-SRP
#  03 - (k,l)-SRP
#  04 - LCSB
# ------------------------------------------------------------------------------
# ./lcsb -alg 1 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
#   -ts ${dPath}.mcs -of ${oFolder}

for m in 64 128 256 512 1024 2048 4096
do
  ./lcsb -alg 2 -n ${n} -qn ${qn} -d ${d} -m ${m} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mcs -of ${oFolder}
done

for m in 64 128 256 512 1024 2048 4096
do
  ./lcsb -alg 4 -n ${n} -qn ${qn} -d ${d} -m ${m} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mcs -of ${oFolder}
done

for k in 10 30 50
do
  for l in 100 200 300
  do
    ./lcsb -alg 3 -n ${n} -qn ${qn} -d ${d} -k ${k} -l ${l} -ds ${dPath}.ds \
      -qs ${dPath}.q -ts ${dPath}.mcs -of ${oFolder}
  done
done