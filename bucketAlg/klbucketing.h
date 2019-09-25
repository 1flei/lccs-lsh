#pragma once

#include "../util.h"
#include <unordered_map>

//bucketing algorithm using (k,L)-setting
//i.e., verify as per one match, from this point of view, this is like L_-inf bucketing
class KLBucketing {
public:
    KLBucketing(int n, int L, int K)
        : nPnts(n)
        , L(L)
        , buckets(L)
        , K(K)
        , cm(nPnts){};
    int nPnts;
    int L;
    std::vector<std::unordered_map<uint32_t, std::vector<int>>> buckets;
    int K;
    CountMarker cm;
    // const std::vector<std::vector<SigType> > *codesp;
    SigType** codesp;

    // std::vector<std::unordered_map<uint64_t, std::vector<int>>> buckets;

    // void build(const std::vector<std::vector<SigType> > &codes) {
    void build(NDArray<2, SigType>& codes)
    {
        assert(codes.lens[1] == K * L);
        // codesp = codes.to_ptr();
        //build hash table based on codes
        for (int i = 0; i < nPnts; i++) {
            for (int j = 0; j < L; j++) {
                uint64_t hashcode = 0;
                for (int k = 0; k < K; k++) {
                    int didx = j * K + k;
                    hashcode = hash_combine(hashcode, uint64_t(codes[i][didx]));
                }

                uint32_t hashcode32 = (hashcode>>7);
                buckets[j][hashcode32].push_back(i);
            }
        }
    }

    void insert(int i, std::vector<SigType>& codei)
    {
        for (int j = 0; j < L; j++) {
            uint64_t hashcode = 0;
            for (int k = 0; k < K; k++) {
                int didx = j * K + k;
                hashcode = hash_combine(hashcode, uint64_t(codei[didx]));
            }
            uint32_t hashcode32 = (hashcode>>7);
            buckets[j][hashcode32].push_back(i);
        }
    }

    template <typename F>
    void for_candidates(int nCandidates, const std::vector<SigType>& qcode, const F& f)
    {
        cm.clear();
        for (int j = 0; j < L; j++) {
            uint64_t hashcode = 0;
            for (int k = 0; k < K; k++) {
                int didx = j * K + k;
                hashcode = hash_combine(hashcode, uint64_t(qcode[didx]));
            }
            uint32_t hashcode32 = (hashcode>>7);
            auto it = buckets[j].find(hashcode32);
            if (it != buckets[j].end()) {
                for (int idx : it->second) {
                    if (!cm.isMarked(idx)) {
                        f(idx);
                        cm.mark(idx);
                        if (!--nCandidates) {
                            return;
                        }
                    }
                }
            }
        }
    }

    int64_t get_memory_usage()
    {
        int64_t ret = sizeof(*this);      //
        ret += sizeof(buckets[0]) * buckets.size(); //the hash_table
        ret += sizeof(uint32_t) * int64_t(nPnts);        //cm && cnt && cntidx

        // printf("ret=%lu\n", ret);

        //hash table
        int64_t entrySize = sizeof(uint32_t) + sizeof(std::vector<int>) + sizeof(void*);
        int64_t bucketSize = sizeof(void*);
        for(auto& ht:buckets){  
            ret += ht.size() * entrySize + ht.bucket_count() * bucketSize;
            // printf("  ret=%lu\n", ret);
            for(auto& p:ht){
                auto& vec = p.second;
                ret += sizeof(int) * vec.capacity();
                // printf("   ret=%lu\n", ret);
            }
        }
        return ret;
    }

};

//bucketing algorithm using (k,L)-setting
//i.e., verify as per one match, from this point of view, this is like L_-inf bucketing
//this version of implementation uses simply vector as hash_table (ignoring the conflication)
class KLBucketingSimpleHasher {
public:
    KLBucketingSimpleHasher(int n, int L, int K)
        : nPnts(n)
        , L(L)
        , buckets(L)
        , K(K)
        , cm(nPnts)
        , hc(K, (uint64_t)nPnts*2)
    {
        //using pertubation hash to combine hash signatures
        for(int i=0;i<L;i++){
            buckets[i].resize(hc.mask+1);
        }
    };
    int nPnts;
    int L;
    std::vector<std::vector<std::vector<int>>> buckets;
    int K;
    CountMarker cm;
    // const std::vector<std::vector<SigType> > *codesp;
    // SigType** codesp;

    HashCombinator<SigType> hc;

    // std::vector<std::unordered_map<uint64_t, std::vector<int>>> buckets;

    // void build(const std::vector<std::vector<SigType> > &codes) {
    void build(NDArray<2, SigType>& codes)
    {
        assert(codes.lens[1] == K * L);
        // codesp = codes.to_ptr();
        //build hash table based on codes
        for (int i = 0; i < nPnts; i++) {
            for (int j = 0; j < L; j++) {
                uint32_t hashcode32 = hc.hash_combine(&codes[i][j*K]);
                buckets[j][hashcode32].push_back(i);
            }
        }
    }

    void insert(int i, std::vector<SigType>& codei)
    {
        for (int j = 0; j < L; j++) {
            uint32_t hashcode32 = hc.hash_combine(&codei[j*K]);
            // assert(hashcode32 < buckets[j].size());
            buckets[j][hashcode32].push_back(i);
        }
    }

    template <typename F>
    void for_candidates(int nCandidates, const std::vector<SigType>& qcode, const F& f)
    {
        cm.clear();
        for (int j = 0; j < L; j++) {
            uint32_t hashcode32 = hc.hash_combine(&qcode[j*K]);
            for (int idx : buckets[j][hashcode32]) {
                if (!cm.isMarked(idx)) {
                    f(idx);
                    cm.mark(idx);
                    if (!--nCandidates) {
                        return;
                    }
                }
            }
        }
    }

    int64_t get_memory_usage()
    {
        int64_t ret = sizeof(*this);      //
        ret += sizeof(buckets[0]) * buckets.size(); //the hash_table
        for(int i=0;i<buckets.size();i++){
            ret += sizeof(buckets[i][0]) * buckets[i].size();
            for(int j=0;j<buckets[i].size();j++){
                ret += sizeof(buckets[i][j]) * buckets[i][j].size();
            }
        }
        ret += sizeof(uint32_t) * int64_t(nPnts);        //cm && cnt && cntidx

        return ret;
    }

};

class KLBucketingCompact {
public:
    KLBucketingCompact(int n, int L, int K)
        : nPnts(n)
        , L(L)
        , buckets(L)
        , K(K)
        , cm(nPnts)
    {
        assert(K<32);
        for(int i=0;i<L;i++){
            buckets[i].resize(1<<K);
        }
    };
    int nPnts;
    int L;
    std::vector<std::vector<std::vector<uint32_t>>> buckets;
    int K;
    CountMarker cm;

    uint64_t get_sig_between(const uint64_t* v, int s, int e)
    {
        assert(e-s<64 && e-s>0);

        uint64_t mask = (1<<(e-s))-1;
        int loc = s/64;
        int shift = s-loc*64;
        return (v[loc]>>shift)&mask;
    }

    // void build(const std::vector<std::vector<SigType> > &codes) {
    void build(NDArray<2, uint64_t>& codes)
    {
        assert(codes.lens[1] == (K*L+63)/64);
        // codesp = codes.to_ptr();
        //build hash table based on codes
        for (int i = 0; i < nPnts; i++) {
            for (int j = 0; j < L; j++) {
                uint32_t codeij = get_sig_between(codes[i], j*K, (j+1)*K);
                buckets[j][codeij].push_back(i);
            }
        }
    }

    void insert(int i, std::vector<uint64_t>& codei)
    {
        for (int j = 0; j < L; j++) {
            uint32_t codeij = get_sig_between(&codei[0], j*K, (j+1)*K);
            buckets[j][codeij].push_back(i);
        }
    }

    template <typename F>
    void for_candidates(int nCandidates, const std::vector<uint64_t>& qcode, const F& f)
    {
        cm.clear();
        for (int j = 0; j < L; j++) {
            uint32_t codeij = get_sig_between(&qcode[0], j*K, (j+1)*K);
            for (int idx : buckets[j][codeij]) {
                if (!cm.isMarked(idx)) {
                    f(idx);
                    cm.mark(idx);
                    if (!--nCandidates) {
                        return;
                    }
                }
            }
        }
    }

    int64_t get_memory_usage()
    {
        int64_t ret = sizeof(*this);      //
        ret += sizeof(buckets[0]) * buckets.size(); //the hash_table
        for(int i=0;i<buckets.size();i++){
            ret += sizeof(buckets[i][0]) * buckets[i].size();
            for(int j=0;j<buckets[i].size();j++){
                ret += sizeof(buckets[i][j]) * buckets[i][j].size();
            }
        }
        ret += sizeof(uint32_t) * int64_t(nPnts);        //cm && cnt && cntidx

        return ret;
    }

};