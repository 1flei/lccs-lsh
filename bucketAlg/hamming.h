#pragma once

#include "../def.h"
#include "../util.h"
#include <vector>
#include "../mih/include/myMihasher.h"

//include linear and MIH

class HammingLinearScan
{
public:
    HammingLinearScan(int d) : dim(d/8) {
        assert(d%8==0);
    };
    int dim;
    int nPnts;
    const uint64_t** codesp;

    std::vector<int> idx;

    void build(NDArray<2, uint64_t> &codes) {
        codesp = codes.to_cptr();
        nPnts  = codes.lens[0];

        // printf("lens=%d, %d,   dim=%d\n", codes.lens[0], codes.lens[1], dim);
        //do nothing but remembering codes & init idx
        idx.resize(codes.lens[0]);
        for(int i=0;i<codes.lens[0];i++){
            idx[i] = i;
        }
    }
    
    template<typename F>
    void for_candidates(int nCandidates, const std::vector<uint64_t> &qcode, const F& f) {
        std::vector<int> dists(nPnts);
        for(int i=0;i<nPnts;i++){
            // dists[i] = calc_ha
            dists[i] = calc_hamming_dist(dim, (uint8_t*)&codesp[i][0], (uint8_t*)&qcode[0]);
        }
        
        const auto& cmp = [&](int a, int b){
            return dists[a] < dists[b];
        };
        nth_element(idx.begin(), idx.begin()+nCandidates, idx.end(), cmp);
        for(int i=0;i<nCandidates;i++){
            printf("hamming-dist[%d]=%d\n", idx[i], dists[idx[i]]);
            f(idx[i]);
        }
    }
    int64_t get_memory_usage()
    {
        return sizeof(*this);
    }
};


//MIH
class MIH
{
public:
    MIH(int d, int m) : dim(d), nChuncks(m) {
        assert(d*sizeof(SigType)*8 / m < 32);
        mih = std::make_unique<myMihasher>(d*sizeof(SigType)*8, m);
        numres.resize(d*sizeof(SigType)*8 + 1);
    };
    int dim;
    int nChuncks;
    // const std::vector<std::vector<SigType> > *codesp;
    std::vector<UINT8*> codesp;

    std::vector<int> idx;

    void build(const std::vector<std::vector<SigType> > &codes) {
        //do nothing but remembering codes & init idx
        codesp.resize(codes.size());
        for(int i=0;i<codesp.size();i++){
            codesp[i] = (UINT8*) &codes[i][0];
        }

        mih->populate(&codesp[0], codes.size());
    }
    
    template<typename F>
    void for_candidates(int nCandidates, const std::vector<SigType> &qcode, const F& f) {
        std::vector<UINT32> results(nCandidates);
        UINT8* u8qcodep = (UINT8*)&qcode[0];
        mih->query(nCandidates, &results[0], &numres[0], u8qcodep);
        for(int i=0;i<nCandidates;i++){
            f(results[i]-1);
        }
    }
protected:
    //as mihasher only accepts flatten binary codes, which does not fit in our framework currently
    //another version using UINT8** is implmented following basically exactly the same structure. (myMihasher)
    std::unique_ptr<myMihasher> mih;
    std::vector<UINT32> numres;
};