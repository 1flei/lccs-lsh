#pragma once

#include <vector>
#include "../util.h"
#include "../pri_queue.h"
#include <random>
#include "../register.h"
#include "../benchmark_util.h"

#include "../hashAlg/srp.h"
#include "../hashAlg/pivots.h"
#include "../bucketAlg/hamming.h"
#include "../bucketAlg/lcs.h"
#include "../bucketAlg/klbucketing.h"

extern bool SRPP_LINEAR_SCAN_REGISTED;
extern bool SRP_LINEAR_SCAN_REGISTED;
extern bool PIVOTSBINARY_LINEAR_SCAN_REGISTED;


//for illustration only
class HasherInterface
{
public:
    // HasherInterface(int d, ...);   //constructor given whatever parameters
    virtual std::vector<SigType> getSig(const Scalar *data);
};
//for illustration only
class BucketerInterface
{
public:
    // BucketerInterface(int k, ...);   //constructor given whatever parameters
    virtual void build(const std::vector<std::vector<SigType> > &codes);

    template<typename F>
    void forCandidates(int nCandidates, const std::vector<SigType> &qcode, const F& f);
};

//simple wraper combining hasher and bucketer
template<typename Hasher, typename Bucketer>
class ComposibleIndex
{
public:
    ComposibleIndex() {}

    //init hasher given whatever parameters
    template<typename... Args>
    void initHasher(Args&&... args){
        hasher = make_unique<Hasher>(std::forward<Args>(args)...);
    }

    template<typename... Args>
    void initBucketer(Args&&... args){
        bucketer = make_unique<Bucketer>(std::forward<Args>(args)...);
    }

    std::unique_ptr<Hasher> hasher;
    std::unique_ptr<Bucketer> bucketer; 

    //use Scalar**
    void build(int n, const Scalar** data){
        assert(hasher && bucketer);

        codes.reserve(n);
        for (int i = 0; i < n; ++i) {
            codes.emplace_back(hasher->getSig(data[i]));
        }
        bucketer->build(codes);
    }

    template<typename FCandidate>
    void query(int nCandidates, const float *query, const FCandidate& f){
        auto qcode = hasher->getSig(query);
        bucketer->forCandidates(nCandidates, qcode, f);
    }

    std::vector<std::vector<SigType> > codes;
};