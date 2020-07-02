#pragma once

#include "../def.h"
#include "../util.h"
#include <vector>

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

