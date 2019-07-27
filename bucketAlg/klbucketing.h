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
        codesp = codes.to_ptr();
        //build hash table based on codes
        for (int i = 0; i < nPnts; i++) {
            for (int j = 0; j < L; j++) {
                uint64_t hashcode = 0;
                for (int k = 0; k < K; k++) {
                    int didx = j * K + k;
                    hashcode = hash_combine(hashcode, uint64_t(codesp[i][didx]));
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
    void forCandidates(int nCandidates, const std::vector<SigType>& qcode, const F& f)
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