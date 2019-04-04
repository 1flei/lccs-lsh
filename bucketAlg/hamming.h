#pragma once

#include "../def.h"
#include "../util.h"
#include <vector>
#include "../mih/include/myMihasher.h"

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


//MIH
class MIH
{
public:
    MIH(int d, int m) : dim(d), nChuncks(m) {
        assert(d*sizeof(SigType)*8 / m < 32);
        mih = make_unique<myMihasher>(d*sizeof(SigType)*8, m);
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
    void forCandidates(int nCandidates, const std::vector<SigType> &qcode, const F& f) {
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