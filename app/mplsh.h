#pragma once

//interface to mplsh (lshkit)
#include "../util.h"
#include "../myndarray.h"
#include <unordered_map>
#include "../lshkit-0.2.1/include/lshkit.h"
// #include <boost/progress.hpp>
#include "../register.h"
#include "../benchmark_util.h"
#include "../hashAlg/e2.h"
#include <queue>

class MPLSH_LSHKIT
{
public:

    typedef lshkit::MultiProbeLshIndex<unsigned> Index;
    typedef std::mt19937 DefaultRng;

    //n: #points; L: #hasher; T: #probes; W: bucket width
    MPLSH_LSHKIT(int n, int dim, int K, int L, float W) 
        : nPnts(n), dim(dim), K(K), L(L), W(W), cm(n)
    {
        Index::Parameter param;

        param.W = W;
        param.range = 1017881; // default value.
        param.repeat = K;
        param.dim = dim;
        DefaultRng rng;

        index.init(param, rng, L);
    };

    int nPnts;
    int dim;
    int K;
    int L;
    float W;
    CountMarker cm;

    const float** data;

    NDArray<2, SigType> codes;
    
    // std::vector<std::unordered_map<uint32_t, std::vector<int>>> buckets;
    Index index;

    // void build(const std::vector<std::vector<SigType> > &codes) {
    void build(const float **data_) 
    {
        data = data_;
        for (unsigned i = 0; i < nPnts; ++i)
        {
            // Insert an item to the hash table.
            // Note that only the key is passed in here.
            // MPLSH will get the feature from the accessor.
            index.insert(i, data[i]);
        }
    }

    void query(unsigned nCandidates, const float* query, MinK_List *list)
    {
        cm.clear();
        // int checkdCnt=0;
        auto scanner = [&](int idx){
            if(!cm.isMarked(idx)){
                cm.mark(idx);
                double dist = calc_l2_dist(dim, data[idx], query);
                list->insert(dist, idx + 1);
            }
        };
        index.query(query, nCandidates, scanner);
    }

    int64_t get_memory_usage()
    {
        //not implemented yet
        return 0;
    }

};