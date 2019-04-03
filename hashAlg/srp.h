#pragma once
#include <vector>
#include "../util.h"
#include <cassert>
#include <random>
#include "../register.h"

//Random Sign Projection
class SRP
{
public:
    SRP(int d, int K, int M);      //dim of data object, #hasher, #bits per hasher 
    ~SRP();
    std::vector<SigType> getSig(const Scalar *data);
protected:
    int dim, K, M;
    std::vector<Scalar> p;
};

//produce K signatures using ~K^0.5 hasher
class SRPPair
{
public:
    SRPPair(int d, int K, int M);      //dim of data object, #hasher, #bits per hasher 
    ~SRPPair();
    std::vector<SigType> getSig(const Scalar *data);
protected:
    int dim, K, M, nHasher;
    std::vector<Scalar> p;
    std::vector<std::pair<int, int> > pairs;
};
