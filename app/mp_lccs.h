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
    MP_LCCS_L2(int n, int dim, int L, float W, double nProbeTimes=2) 
        : nPnts(n), dim(dim), L(L), W(W), cm(n), 
        lookupLen(std::min(L-1, lookupLen) ), nProbeTimes(nProbeTimes), 
        hasher(dim, L, W, nProbeTimes), bucketer(L, 1)
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
        int nProbes = nProbeTimes*L+1;
        hasher.forSig(nProbes, query, [&](const std::vector<E2MP::SigType>& qcode, int firstIdx){
            if(firstIdx==-1){
                // printf("hq0\n");
                // printVec(qcode.begin(), L);
                bucketer.for_candidates(checkedCandidates, qcode, f);
                // printf("hq0 done!!!\n");
            } else{
                int startIdx = (firstIdx-lookupLen+L+1)%L;
                //let startIdx be the first idx such that matched_loc will be affected
                for(int i=1;i<16;i++){
                    int curloc = (firstIdx-i+L)%L;
                    // printf("   %d, %d, %d\n", i, bucketer.res_lowlen[curloc], bucketer.res_highlen[curloc]);
                    if(bucketer.res_lowlen[curloc] < i || bucketer.res_highlen[curloc] < i){
                        startIdx = curloc;
                        break;
                    }
                }
                // printf("hq[%d, %d]!!!\n", startIdx, firstIdx);
                bucketer.for_candidates_between(startIdx, firstIdx+1, checkedCandidates, qcode, f);
                // printf("hq[%d, %d] DONE!!!\n", startIdx, firstIdx);
            }
        });
    }


    int nPnts;
    int dim;
    int L;
    float W;
    CountMarker cm;
    int lookupLen;
    double nProbeTimes;

    const float** data;

    NDArray<2, SigType> codes;

    E2MP hasher;
    LCCS_SORT_INT bucketer; 
};

//combination of LCCS and E2MP
class MP_LCCS_CP
{
public:
    MP_LCCS_CP(int n, int dim, int L, double nProbeTimes=2) 
        : nPnts(n), dim(dim), L(L), cm(n), 
        lookupLen(std::min(L-1, lookupLen) ), nProbeTimes(nProbeTimes), 
        hasher(dim, L, nProbeTimes), bucketer(L, 1)
    {
    };
    int nPnts;
    int dim;
    int L;
    CountMarker cm;
    int lookupLen;
    double nProbeTimes;
    CrossPolytopeMP hasher;

    const float** data;

    NDArray<2, SigType> codes;

    LCCS_SORT_INT bucketer; 
    
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
        int nProbes = nProbeTimes*L+1;
        hasher.forSig(nProbes, query, [&](const std::vector<E2MP::SigType>& qcode, int firstIdx){
            if(firstIdx==-1){
                // printf("hq0\n");
                // printVec(qcode.begin(), L);
                bucketer.for_candidates(checkedCandidates, qcode, f);
                // printf("hq0 done!!!\n");
            } else{
                int startIdx = (firstIdx-lookupLen+L+1)%L;
                //let startIdx be the first idx such that matched_loc will be affected
                for(int i=1;i<16;i++){
                    int curloc = (firstIdx-i+L)%L;
                    // printf("   %d, %d, %d\n", i, bucketer.res_lowlen[curloc], bucketer.res_highlen[curloc]);
                    if(bucketer.res_lowlen[curloc] < i || bucketer.res_highlen[curloc] < i){
                        startIdx = curloc;
                        break;
                    }
                }
                // printf("hq[%d, %d]!!!\n", startIdx, firstIdx);
                bucketer.for_candidates_between(startIdx, firstIdx+1, checkedCandidates, qcode, f);
                // printf("hq[%d, %d] DONE!!!\n", startIdx, firstIdx);
            }
        });
    }


};

}
