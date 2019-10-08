#pragma once
#include "../util.h"
#include "../register.h"
#include "../benchmark_util.h"
#include <queue>
#include "falconn/lsh_nn_table.h"
#include "falconn/core/polytope_hash.h"

//meaning cross_polytope + MPLSH for cosine similarity
class FALCONN_INDEX
{
public:
    typedef falconn::DenseVector<Scalar> PointType;
    typedef Eigen::Map<PointType> MPoint;

    // std::unique_ptr<falconn::LSHNearestNeighborQuery<PointType>> queryObject;
    const float** data;

    int nPnts;
    int dim;
    int L;

    FALCONN_INDEX(int n, int d, int L, int K, int nRotations=1, bool useBitPackedHashTable=false) 
        : nPnts(n), dim(d), L(L)
    {
        params.dimension = d;
        params.lsh_family = falconn::LSHFamily::CrossPolytope;
        params.l = L;
        params.distance_function = falconn::DistanceFunction::EuclideanSquared;
        params.num_rotations = nRotations;

        if(!useBitPackedHashTable){
            params.k = K;
            params.last_cp_dimension = d;
            params.storage_hash_table = falconn::StorageHashTable::LinearProbingHashTable;
        } else {
            //if use bitPackedHash, let k be the num of bits used instead of the number of hashers cancatenated
            //and let FALCONN itself setting the number of hashers and last_cp_dim
            int nBits = K;
            falconn::compute_number_of_hash_functions<PointType>(nBits, &params);

            printf("calculated_k=%d, last_cp_dim=%d\n", params.k, params.last_cp_dimension);
            params.storage_hash_table = falconn::StorageHashTable::BitPackedFlatHashTable;
        }

        // falconn::compute_number_of_hash_functions<PointType>(nBits, &params);
        // printf("num_hash_functions=%d, last_cp_dimensions=%d\n", params.k, params.last_cp_dimension);
        params.num_setup_threads = 1;
    }
    
    std::vector<PointType> data_vec;

    void build(const float **data_) 
    {
        data = data_;
        for(int i=0;i<nPnts;i++){
            //will copy
            data_vec.push_back(MPoint((float*)data[i], dim));
        }
        table = falconn::construct_table<PointType>(data_vec, params);
    }
    void query(unsigned nProbes, const float* query, MinK_List *list)
    {
        int maxNumberCandidates = nPnts/10;
        std::unique_ptr<falconn::LSHNearestNeighborQuery<PointType>> queryObjectPtr =
            table->construct_query_object(nProbes+L, maxNumberCandidates);

        std::vector<int32_t> candidates;
        queryObjectPtr->get_unique_candidates(MPoint((float*)query, dim), &candidates);
        for(int idx:candidates){
            float l2dist = calc_angle(dim, data[idx], query);
            // printf("%d, %f\n", idx, l2dist);
            list->insert(l2dist, idx + 1);
        }
    }

    int64_t get_memory_usage()
    {
        return 0;
    }

    falconn::LSHConstructionParameters params;
    std::unique_ptr<falconn::LSHNearestNeighborTable<PointType> > table;
};