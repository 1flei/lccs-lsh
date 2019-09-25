#pragma once

#include "../benchmark_util.h"
#include "../myndarray.h"
#include "../pri_queue.h"
#include "../register.h"
#include "../util.h"
#include <random>
#include <vector>

#include "../bucketAlg/c2bucketing.h"
#include "../bucketAlg/hamming.h"
#include "../bucketAlg/klbucketing.h"
#include "../bucketAlg/lcs.h"
#include "../bucketAlg/lcs_int.h"
#include "../bucketAlg/lcs_int_reorder.h"
#include "../bucketAlg/mymplsh.h"
// #include "../bucketAlg/lcs_opt.h"
// #include "../bucketAlg/lcs_test.h"
#include "../hashAlg/e2.h"
#include "../hashAlg/pivots.h"
#include "../hashAlg/srp.h"
#include "../hashAlg/e2eigen.h"
#include "../hashAlg/polytope.h"

//for illustration only
class HasherInterface {
public:
    // HasherInterface(int d, ...);   //constructor given whatever parameters
    virtual std::vector<SigType> getSig(const Scalar* data);
    virtual void getSig(const Scalar* data, SigType* sig);
};
//for illustration only
class BucketerInterface {
public:
    // BucketerInterface(int k, ...);   //constructor given whatever parameters
    virtual void build(const std::vector<std::vector<SigType>>& codes);

    template <typename F>
    void for_candidates(int nCandidates, const std::vector<SigType>& qcode, const F& f);
};

//simple wraper combining hasher and bucketer
//assume continuous data arrangement
//e.g. both data, query and hash signature will be represented as 2darray
template <typename Hasher, typename Bucketer, typename SigT = SigType>
class ComposibleIndex {
public:
    ComposibleIndex() {
        codes.resize({size_t(0), size_t(0)});
    }

    //init hasher given whatever parameters
    template <typename... Args>
    void initHasher(Args&&... args)
    {
        hasher = std::make_unique<Hasher>(std::forward<Args>(args)...);
    }

    template <typename... Args>
    void initBucketer(Args&&... args)
    {
        bucketer = std::make_unique<Bucketer>(std::forward<Args>(args)...);
    }

    std::unique_ptr<Hasher> hasher;
    std::unique_ptr<Bucketer> bucketer;

    const double printStep = 0.1;

    //use Scalar**
    void build(int n, const Scalar** data)
    {
        assert(hasher && bucketer);

        // codes.reserve(n);
        printf("computing hash-code\n");
        codes.resize({ size_t(n), size_t(hasher->sigdim) });
        for (int i = 0; i < n; ++i) {
            // codes.emplace_back(hasher->getSig(data[i]));
            hasher->getSig(data[i], &codes[i][0]);
        }
        printf("building index\n");
        bucketer->build(codes);
    }

    //build index and then free sigs
    void build_wo_sigs(int n, const Scalar** data)
    {
        assert(hasher && bucketer);

        // codes.reserve(n);
        printf("computing hash-code\n");
        codes.resize({ size_t(n), size_t(hasher->sigdim) });
        for (int i = 0; i < n; ++i) {
            // codes.emplace_back(hasher->getSig(data[i]));
            hasher->getSig(data[i], &codes[i][0]);
        }
        printf("building index\n");
        bucketer->build(codes);

        printf("free codes\n");
        codes.resize({0, 0});
    }

    void inc_build(int n, const Scalar** data)
    {
        assert(hasher && bucketer);

        int printcnt = 0;
        printf("incremental building\n");
        for (int i = 0; i < n; i++) {
            if (i %1000==0) {
                double progress = i*1. / n;
                print_progress(progress);
            }
            auto codei = hasher->getSig(data[i]);
            bucketer->insert(i, codei);
        }
        print_progress(1.);
        printf("\n");
    }

    void print_progress(double percentage)
    {
        const char PBSTR[] = "====================";
        const int PBWIDTH = 20;
        int val = (int)(percentage * 100);
        int lpad = (int)(percentage * PBWIDTH);
        int rpad = PBWIDTH - lpad;
        printf("\r%3d%% [%.*s>%*s]", val, lpad, PBSTR, rpad, "");
        fflush(stdout);
    }

    int64_t get_memory_usage()
    {
        assert(hasher && bucketer);
        //code
        int64_t codes_usage = codes.get_memory_usage();
        printf("codes_usage=%ld\n", codes_usage);
        int64_t hash_usage = hasher->get_memory_usage();
        printf("hash_usage=%ld\n", hash_usage);
        //bucket index
        int64_t bucket_usage = bucketer->get_memory_usage();
        printf("bucket_usage=%ld\n", bucket_usage);

        return codes_usage + hash_usage + bucket_usage;
    }

    template <typename FCandidate>
    void query(int nCandidates, const float* query, const FCandidate& f)
    {
        auto qcode = hasher->getSig(query);
        bucketer->for_candidates(nCandidates, qcode, f);
    }

    // std::vector<std::vector<SigT> > codes;
    NDArray<2, SigT> codes;
};