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
    MP_LCCS_L2(int n, int dim, int L, float W, int lookupLen=8, int checkedCandidates=4) 
        : nPnts(n), dim(dim), L(L), W(W), cm(n), 
        lookupLen(lookupLen), checkedCandidates(checkedCandidates), 
        hasher(dim, L, W), bucketer(L, 1)
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
    void query(int nProbes, const float* query, const FCandidate& f)
    {
        hasher.forSig(nProbes, query, [&](const std::vector<E2MP::SigType>& qcode, int lastIdx){
            int startIdx = (lastIdx-lookupLen+L)%L;
            bucketer.for_candidates_between(startIdx, lastIdx, checkedCandidates, qcode, f);
        });
    }


    int nPnts;
    int dim;
    int L;
    float W;
    CountMarker cm;
    int lookupLen;
    int checkedCandidates;

    const float** data;

    NDArray<2, SigType> codes;

    E2MP hasher;
    LCCS_SORT_INT bucketer; 
};

}
