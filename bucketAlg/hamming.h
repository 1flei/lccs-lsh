#pragma once

#include "../def.h"
#include "../util.h"
#include <vector>

//include linear and MIH



class HammingLinearScan
{
public:
    HammingLinearScan(int d) : dim(d) {};
    int dim;
    const std::vector<std::vector<SigType> > *codesp;

    std::vector<int> idx;

    void build(const std::vector<std::vector<SigType> > &codes) {
        //do nothing but remembering codes & init idx
        codesp = &codes;

        idx.resize(codes.size());
        for(int i=0;i<codes.size();i++){
            idx[i] = i;
        }
    }
    
    template<typename F>
    void forCandidates(int nCandidates, const std::vector<SigType> &qcode, const F& f) {
        const auto& codes = *codesp;
        std::vector<int> dists(codes.size());
        for(int i=0;i<codes.size();i++){
            // dists[i] = calc_ha
            dists[i] = calc_hamming_dist(dim * sizeof(SigType), (const uint8_t*) &codes[i][0], (const uint8_t*) &qcode[0]);
        }
        
        const auto& cmp = [&](int a, int b){
            return dists[a] < dists[b];
        };
        nth_element(idx.begin(), idx.begin()+nCandidates, idx.end(), cmp);
        for(int i=0;i<nCandidates;i++){
            f(idx[i]);
        }
    }
};