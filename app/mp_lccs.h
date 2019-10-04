#pragma once

#include "../hashAlg/mp_hasher.h"
#include "../myndarray.h"
#include "../bucketAlg/lcs_int.h"
#include "../register.h"
#include "../benchmark_util.h"

namespace mylccs
{

//combination of LCCS and E2MP
class MP_LCCS_L2
{
public:
    MP_LCCS_L2(int n, int dim, int L, float W, int lookupLen=7, int defaultNProbe=2) 
        : nPnts(n), dim(dim), L(L), W(W), cm(n), 
        lookupLen(std::min(L-1, lookupLen) ), defaultNProbe(defaultNProbe), 
        hasher(dim, L, W, lookupLen), bucketer(L, 1)
    {
    };
    
    void build(const Scalar** data)
    {
        // codes.reserve(n);
        printf("computing hash-code\n");
        codes.resize({ size_t(nPnts), size_t(hasher.sigdim) });
        for (int i = 0; i < nPnts; ++i) {
            // codes.emplace_back(hasher->getSig(data[i]));
            hasher.getSig(data[i], &codes[i][0]);
        }
        printf("building index\n");
        bucketer.build(codes);
    }

    int64_t get_memory_usage()
    {
        //not implemented yet
        return 0;
    }

    template <typename FCandidate>
    void query(int checkedCandidates, const float* query, const FCandidate& f)
    {
        // int checkedCandidates = 4*nProbes;
        // int checkedCandidates = 16;
        int nProbes = defaultNProbe*L+1;
        hasher.forSig(nProbes, query, [&](const std::vector<E2MP::SigType>& qcode, int lastIdx){
            if(lastIdx==-1){
                // printf("hq0\n");
                // printVec(qcode.begin(), L);
                bucketer.for_candidates(checkedCandidates, qcode, f);
                // printf("hq0 done!!!\n");
            } else{
                int startIdx = (lastIdx-lookupLen+L)%L;
                // printf("hq[%d, %d]\n", startIdx, lastIdx);
                // printVec(qcode.begin(), L);
                bucketer.for_candidates_between(startIdx, lastIdx, checkedCandidates, qcode, f);
                // printf("hq[%d, %d] done!!!\n", startIdx, lastIdx);
            }
        });
    }


    int nPnts;
    int dim;
    int L;
    float W;
    CountMarker cm;
    int lookupLen;
    int defaultNProbe;

    const float** data;

    NDArray<2, SigType> codes;

    E2MP hasher;
    LCCS_SORT_INT bucketer; 
};

}
