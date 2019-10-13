#pragma once

#include <vector>
#include "../util.h"
#include <cassert>
#include <random>


//Random Projection
class E2
{
public:
    using SigType = int32_t;
    E2(int d, int K, double r)      //dim of data object, #hasher, radius 
        :dim(d), K(K), r(r)
    {
        assert(d > 0 && K > 0);

        std::normal_distribution<double> normal(0.);
        std::uniform_real_distribution<double> uniform(0., r);
        // std::random_device rd;
        // std::default_random_engine rng(rd());
        std::default_random_engine rng(GLOBAL_SEED);

        p.resize(K*d);
        b.resize(K);
        for (int i = 0; i < K * d; i++) {
            p[i] = normal(rng);
        }
        for (int i = 0; i < K; i++) {
            b[i] = uniform(rng);
        }
        sigdim = K;
    }
    ~E2() {}
    std::vector<SigType> getSig(const Scalar *data) const
    {
        std::vector<SigType> ret(sigdim);
        getSig(data, &ret[0]);
        return ret;
    }
    void getSig(const Scalar *data, SigType* ret) const
    {
        for(int k=0;k<K;k++){
            double projection = 0.;
            for(int i=0;i<dim;i++){
                double x = data[i];
                int pidx = k*dim + i;
                projection += x*p[pidx];
            }
            projection += b[k];

            ret[k] = SigType(floor(projection/r) );
        }
    }

    void getSigDelta(const Scalar *data, SigType* ret, Scalar* delta) const
    {
        for(int k=0;k<K;k++){
            double projection = 0.;
            for(int i=0;i<dim;i++){
                double x = data[i];
                int pidx = k*dim + i;
                projection += x*p[pidx];
            }
            projection += b[k];

            ret[k] = SigType(floor(projection/r) );

            //only for MPLSH
            delta[k] = projection/r - ret[k];
        }
    }

    int64_t get_memory_usage()
    {
        return int64_t(sizeof(*this)) + int64_t(p.capacity())*sizeof(Scalar) + int64_t(b.capacity())*sizeof(Scalar);
    }

    int dim, K;
    double r;
    int sigdim;
protected:
    std::vector<Scalar> p;
    std::vector<Scalar> b;
};