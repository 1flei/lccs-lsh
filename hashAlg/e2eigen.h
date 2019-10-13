#pragma once

#include <vector>
#include "../util.h"
#include <cassert>
#include <random>
#include "../Eigen/Eigen"


//Random Projection using eigen
class E2Eigen
{
public:
    using SigType = int32_t;
    E2Eigen(int d, int K, double r, int nComb=1)      //dim of data object, #hasher, radius 
        :dim(d), K(K), r(r)
    {
        assert(d > 0 && K > 0);

        std::normal_distribution<double> normal(0.);
        std::uniform_real_distribution<double> uniform(0., r);
        // std::random_device rd;
        // std::default_random_engine rng(rd());
        std::default_random_engine rng(GLOBAL_SEED);

        p.resize(K, d);
        b.resize(K);
        for (int i = 0; i < K; i++) {
            for(int j=0;j<d;j++){
                p(i, j) = normal(rng);
            }
        }
        for (int i = 0; i < K; i++) {
            b(i) = uniform(rng) + 0.5*r;
        }
        sigdim = K;
    }
    ~E2Eigen() {}
    std::vector<SigType> getSig(const Scalar *data)
    {
        std::vector<SigType> ret(sigdim);
        getSig(data, &ret[0]);
        return ret;
    }
    void getSig(const Scalar *data, SigType* ret)
    {
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);
        Eigen::Map<Eigen::Matrix<SigType, Eigen::Dynamic, 1> > ret_v(ret, sigdim);

        auto ret_r = (p*data_v+b)/r;
        ret_v = ret_r.cast<int32_t>();
    }

    int64_t get_memory_usage()
    {
        return int64_t(sizeof(*this)) + (K*dim)*sizeof(Scalar) + int64_t(K)*sizeof(Scalar);
    }

    int dim, K;
    Scalar r;
    int sigdim;
protected:
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> p;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b;
    
    // std::vector<Scalar> p;
    // std::vector<Scalar> b;
};