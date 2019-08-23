#pragma once

#include "../util.h"
#include <unordered_map>
#include "../myndarray.h"

class StepWiseProbeGen
{
public:
    template<typename VecIT>
    StepWiseProbeGen(const VecIT& beg, const VecIT& end) : code(beg, end)
    {
    }

    std::vector<SigType> code;
    std::vector<uint32_t> code_hashs;
    std::hash<SigType> hasher;

    template<typename F> 
    bool stepwise_r(uint64_t hashcode, int curloc, int dist, const F& f)
    {
        // printf("%d, %d\n", curloc, dist);
        if(curloc >= code.size()){
            uint32_t hashcode32 = (hashcode>>7);
            return f(hashcode32);
        }
        //not change
        uint64_t hashcode_0 = hash_combine(hashcode, uint64_t(code[curloc]));
        if(stepwise_r(hashcode_0, curloc+1, dist, f)){
            return true;
        }

        for(int i=1;i<=dist;i++){
            uint64_t hashcode_ip = hash_combine(hashcode, uint64_t(code[curloc]+i));
            if(stepwise_r(hashcode_ip, curloc+1, dist-i, f) ){
                return true;
            }

            uint64_t hashcode_in = hash_combine(hashcode, uint64_t(code[curloc]-i));
            if(stepwise_r(hashcode_in, curloc+1, dist-i, f) ){
                return true;
            }
        }

        return false;
    }

    //f :: uint32_t -> bool
    template<typename F>
    void forProbes(const F& f){
        for(int d=0;;d++){
            // printf("d=%d\n", d);
            if(stepwise_r(0, 0, d, f)){
                return ;
            }
        }
    }
};


class MYMPLSH
{
public:
    //n: #points; L: #hasher; T: #probes; W: bucket width
    MYMPLSH(int n, int L, int K)
        : nPnts(n)
        , K(K)
        , L(L)
        , buckets(L)
        , cm(nPnts){};
    int nPnts;
    int K;
    int L;
    std::vector<std::unordered_map<uint32_t, std::vector<int>>> buckets;
    CountMarker cm;
    // const std::vector<std::vector<SigType> > *codesp;
    SigType** codesp;

    typedef StepWiseProbeGen ProbeGen;
    std::hash<SigType> hasher;



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
    void for_candidates(int nProbes, const std::vector<SigType>& qcode, const F& f)
    {
        cm.clear();
        for (int j = 0; j < L; j++) {
            int nProbes_j = nProbes;
            ProbeGen swpg(&qcode[j*K], &qcode[(j+1)*K]);
            swpg.forProbes([&](uint32_t hashcode32) -> bool{
                auto it = buckets[j].find(hashcode32);
                if (it != buckets[j].end()) {
                    // printf("nCandidates=%d\n", it->second.size());
                    for (int idx : it->second) {
                        if (!cm.isMarked(idx)) {
                            f(idx);
                            cm.mark(idx);
                        }
                    }
                }
                if(--nProbes_j<=0){
                    return true;
                } else{
                    return false;
                }
            });
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