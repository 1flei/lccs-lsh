#pragma once
#include <vector>
#include "../util.h"
#include <cassert>
#include <random>
#include "../register.h"
#include "../FFHT/fht_header_only.h"
//Fast Hadamard Transform implementation

//Apply Hadamard transform 3 times 
//x |-> H D_3 H D_2 H D_1 x produces d signatures             D_i is a random disgonal +-1-matrix for all i
//cost O(d) space and O(dlogd) time using Fast Hadamard Transform to produce d bits signatures
//thus to produce K signatures with M bits, run K*M/d independent one
//in total use O(KM) space and O(KMlog(d)) time
//refers to [KDD'11] "Fast Locality-Sensitive Hashing" & [NIPS'15] "Practical and Optimal LSH for Angular Distance" for more detail
class HD3Hash
{
public:
    HD3Hash(int d, int K, int M);      //dim of data object, #hasher, #bits per hasher 
    ~HD3Hash();
    std::vector<SigType> getSig(const Scalar *data);
protected:
    int dim, K, M, nHasher;
    std::vector<Scalar> p;
    std::vector<std::pair<int, int> > pairs;
};