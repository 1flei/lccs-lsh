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
    void getSig(const Scalar *data, SigType* sig);

    int dim, K, M;
    int sigdim;
protected:
    std::vector<Scalar> p;
};

//Random Sign Projection
class SRPCompact
{
public:
    SRPCompact(int d, int K);      //dim of data object, #hasher
    ~SRPCompact();
    std::vector<uint64_t> getSig(const Scalar *data);
    void getSig(const Scalar *data, uint64_t* sig);

    int dim, K;
    int sigdim;

    int64_t get_memory_usage()
    {
        return int64_t(sizeof(*this)) + sizeof(p) + sizeof(Scalar)*p.capacity();
    }

protected:
    std::vector<Scalar> p;
};

