#pragma once

#include "../util.h"
#include <unordered_map>

//bucketing algorithm using (k,L)-setting
//i.e., verify as per one match, from this point of view, this is like L_-inf bucketing
class C2Bucketing
{
public:
    // C2Bucketing(int n, int L, int K) : nPnts(n), L(L), buckets(L), K(K), cm(0){};
    C2Bucketing(int n, int L) : nPnts(n), L(L), buckets(L), cm(0){
        cm.resize(nPnts);
        cnt.resize(nPnts);
        cntidx.resize(nPnts);
        for(int i=0;i<nPnts;i++){
            cntidx[i] = i;
        }
    };
    int nPnts;
    int L;
    int K;      //threshold
    // const std::vector<std::vector<SigType> > *codesp;
    SigType** codesp;

    std::vector<std::unordered_map<SigType, std::vector<int> > > buckets;
    CountMarker cm;
    std::vector<int> cnt;
    std::vector<int> cntidx;

    // void build(const std::vector<std::vector<SigType> > &codes) {
    void build(NDArray<2, SigType> &codes) {
        assert(codes.lens[1] == L);

        // buckets.reserve(L);
        codesp = codes.to_ptr();
        //build hash table based on codes
        for(int i=0;i<nPnts;i++){
            for(int j=0;j<L;j++){
                SigType codeij = codesp[i][j];
                buckets[j][codeij].push_back(i);
            }
        }
    }

    void insert(int i, std::vector<SigType>& codei)
    {
        for (int j = 0; j < L; j++) {
            SigType codeij = codei[j];
            buckets[j][codeij].push_back(i);
        }
    }

    //cnt strategy, it turns out that cnt strategy is faster than threshold strategy
    template<typename F>
    void forCandidates(int nCandidates, const std::vector<SigType> &qcode, const F& f) 
    {
        // int nMarked = 0;    
        // cm.clear();
        buckets.reserve(L);
        std::fill(cnt.begin(), cnt.end(), 0);
        for(int j=0;j<L;j++){
            SigType qcodej = qcode[j];
            auto it = buckets[j].find(qcodej);
            if(it !=buckets[j].end()){
                for(int idx:it->second){
                    cnt[idx]++;
                }
            }
        }

        std::sort(cntidx.begin(), cntidx.end(), [&](int a, int b){
            return cnt[a] > cnt[b];
        });
        for(int i=0;i<nCandidates;i++){
            f(cntidx[i]);
        }
    }

    int64_t get_memory_usage()
    {
        int64_t ret = sizeof(*this);      //the vector
        ret += sizeof(buckets[0]) * buckets.size(); //the hash_table
        ret += sizeof(uint32_t) * int64_t(nPnts) + sizeof(int) * int64_t(nPnts) *2;        //cm && cnt && cntidx

        // printf("ret=%ld\n", ret);
        //hash table
        int64_t entrySize = sizeof(SigType) + sizeof(std::vector<int>) + sizeof(void*);
        int64_t bucketSize = sizeof(void*);
        for(auto& ht:buckets){  
            ret += ht.size() * entrySize + ht.bucket_count() * bucketSize;
            // printf("  ret=%ld\n", ret);
            for(auto& p:ht){
                auto& vec = p.second;
                ret += sizeof(int) * vec.capacity();
                // printf("    ret=%ld, cap=%d\n", ret, vec.capacity());
            }
        }
        return ret;
    }

};