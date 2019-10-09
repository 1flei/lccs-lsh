#pragma once

#include "falconn/lsh_nn_table.h"
#include "falconn/core/polytope_hash.h"
#include "../util.h"

//data is assumed to be normalized
//a wrapper of polytope hasher wrapper using falconn library
class PolytopeHasher
{
public:
    // using SigType = uint32_t;
    // using Hasher = falconn::core::CrossPolytopeHashBase<falconn::core::CrossPolytopeHashDense<Scalar, SigType>, std::vector<float>, Scalar, SigType>;
    using Hasher = falconn::core::CrossPolytopeHashDense<Scalar, SigType>;
    using Transformer = Hasher::HashTransformation;
    using TransformedVectorType = Hasher::TransformedVectorType;

    int dim;
    // int l;
    int sigdim;
    int num_rotations;
    int last_cp_dim;
    Hasher hasher;
    Transformer transfromer;

    PolytopeHasher(int vector_dim,
                int l, int num_rotations=1,
                int seed=666)
        :dim(vector_dim), 
         sigdim(l), 
         num_rotations(num_rotations), 
         last_cp_dim(vector_dim), 
         //k=1
         hasher(vector_dim, 1, l, num_rotations, last_cp_dim, seed), 
         transfromer(hasher)
    {
        hasher.reserve_transformed_vector_memory(&transformedVec);
    }

    std::vector<SigType> getSig(const Scalar *data)
    {
        std::vector<SigType> ret(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);

        // Eigen::Matrix<Scalar, Eigen::Dynamic, 1> normalized_data_v = data_v.normalized();
        hasher.hash(data_v, &ret, &transformedVec);
        return ret;
    }

    void getSig(const Scalar *data, SigType* ret)
    {
        std::vector<SigType> res(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);

        // Eigen::Matrix<Scalar, Eigen::Dynamic, 1> normalized_data_v = data_v.normalized();
        hasher.hash(data_v, &res, &transformedVec);
        std::copy(res.begin(), res.end(), ret);
    }

    TransformedVectorType transformedVec;

    int64_t get_memory_usage()
    {
        return 0;
    }
};