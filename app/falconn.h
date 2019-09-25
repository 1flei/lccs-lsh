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

    FALCONN_INDEX(int n, int d, int L, int nBits, int nRotations=1) 
        : nPnts(n), dim(d), L(L)
    {
        params.dimension = d;
        params.lsh_family = falconn::LSHFamily::CrossPolytope;
        params.l = L;
        params.distance_function = falconn::DistanceFunction::EuclideanSquared;
        falconn::compute_number_of_hash_functions<PointType>(nBits, &params);
        printf("num_hash_functions=%d\n", params.k);

        params.num_rotations = nRotations;
        params.num_setup_threads = 1;
        params.storage_hash_table = falconn::StorageHashTable::BitPackedFlatHashTable;
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
            list->insert(l2dist, idx + 1);
        }
    }

    int64_t get_memory_usage()
    {
        //not implemented yet
        //this class
        int64_t ret = sizeof(*this);
        //table, may ignore some insignificant part such as some pointers
        //table is essentially LSHNNTableWrapper pointer
        //memory_usage
        using namespace falconn;
        using KeyType = int32_t;
        typedef typename PointTypeTraits<PointType>::ScalarType ScalarType;

        typedef typename wrapper::PointTypeTraitsInternal<PointType>::CosineDistance
            DistanceFunctionType;

        typedef uint32_t HashType;

        typedef typename wrapper::PointTypeTraitsInternal<
            PointType>::template HPHash<HashType>
            LSHType;
        typename std::unique_ptr<LSHType> LSHPointerType;
        typedef core::BitPackedFlatHashTable<HashType> HashTable;
        typedef typename HashTable::Factory HashTableFactoryType;
        typedef std::vector<PointType> PointSet;

        typedef typename wrapper::DataStorageAdapter<PointSet>::template DataStorage<KeyType>
            DataStorageType;

        typedef core::StaticCompositeHashTable<HashType, KeyType, HashTable>
            CompositeHashTableType;
        typedef core::StaticLSHTable<PointType, KeyType, LSHType, HashType,
                                     CompositeHashTableType, DataStorageType>
            LSHTableType;

        using TableWrapper = falconn::wrapper::LSHNNTableWrapper<PointType, KeyType, ScalarType,
                                    DistanceFunctionType, LSHTableType,
                                    LSHType, HashTableFactoryType,
                                    CompositeHashTableType, DataStorageType>;
        
        //hacker ;)
        // struct MyWrapper : public TableWrapper
        // {
        //     MyWrapper(const TableWrapper& tw):
        //         TableWrapper(tw) {}

        //     int64_t get_memory_usage() {
        //         int64_t ret = 0;
        //         //those are the member that we need to count the memory usage
        //         // lsh_, lsh_table_, hash_table_factory_, composite_hash_table_, data_storage_;
        //         ret += data_storage_->size() * sizeof(KeyType);
        //         return ret;
        //     }
        // };
        
        // MyWrapper* table_tmp = dynamic_cast<MyWrapper*>(table.get());
        // ret += table_tmp->get_memory_usage();
        return 0;
    }

    falconn::LSHConstructionParameters params;
    std::unique_ptr<falconn::LSHNearestNeighborTable<PointType> > table;
};