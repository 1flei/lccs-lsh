#pragma once

//interface to mplsh (lshkit)
#include "../util.h"
#include <unordered_map>
#include "../lshkit-0.2.1/include/lshkit.h"
#include <boost/progress.hpp>
#include "../register.h"
#include "../benchmark_util.h"
#include "../bucketAlg/mymplsh.h"
#include "../hashAlg/e2.h"
#include <queue>

struct Score
{
    int loc;
    float score;
    int direction;
    Score(int loc, int score, int direction): loc(loc), score(score), direction(direction) {}
    Score() {}

    bool operator<(const Score& s) const {
        return score < s.score;
    }
};

struct Probe
{
    Probe(int idx, float score) 
        : mask(1<<idx), tot_score(score), maxidx(idx) 
    {}

    uint64_t mask;
    float tot_score;
    float last_score;
    int maxidx;

    // template<typename IT>
    // uint32_t hash_combine_reduce(const IT& beg, const IT& end)
    // {
    //     uint64_t hashcode = 0;
    //     for (auto it=beg;it!=end;++it) {
    //         hashcode = hash_combine(hashcode, uint64_t(*it));
    //     }
    //     uint32_t hashcode32 = (hashcode>>7);
    //     return hashcode32;
    // }

    bool operator < (const Probe &p) const { return tot_score+last_score < p.tot_score+p.last_score; }
    uint32_t getHash(int beg, int end, const std::vector<SigType>& qcode, const std::vector<Score>& scores)
    {
        std::vector<SigType> qcode_copy(qcode.begin()+beg, qcode.begin()+end);
        for(int i=0;i<scores.size();i++){
            if((1<<i) & mask){
                qcode_copy[scores[i].loc] -= scores[i].direction;
            }
        }

        // printf("mask=%x\n", mask);
        // printVec(qcode_copy.begin(), qcode_copy.end(), "qcode_copy:");

        // hash_combine_reduce(qcode_copy.begin(), qcode_copy.end());
        uint64_t hashcode = 0;
        for (int k = 0; k < qcode_copy.size(); k++) {
            hashcode = hash_combine(hashcode, uint64_t(qcode_copy[k]));
        }
        uint32_t hashcode32 = (hashcode>>7);
        return hashcode32;
    }

    void expand(float new_score)
    {
        tot_score += last_score;
        last_score = new_score;

        ++maxidx;
        mask |= (1<<maxidx);
    }

    void shift(float new_score)
    {
        last_score = new_score;

        mask &= ~(1<<maxidx);
        ++maxidx;
        mask |= (1<<maxidx);
    }
}; 


class FancyProbeGen
{
public:
    //FancyProbeGen may need to know the information of raw data
    FancyProbeGen(int K, int L, const Scalar* query, const E2& e2lsh)
        :K(K), L(L), e2lsh(e2lsh), qcode(e2lsh.sigdim), deltas(e2lsh.sigdim)
    {
        e2lsh.getSigDelta(query, &qcode[0], &deltas[0]);

        scores.resize(L);
        for (int j = 0; j < L; j++) {
            scores[j].resize(2*K);
            for (int k = 0; k < K; k++) {
                int didx = j * K + k;
                scores[j][2*k].loc = k;
                scores[j][2*k].direction = 1;    // direction
                scores[j][2*k].score = deltas[didx];
                scores[j][2*k+1].loc = k;
                scores[j][2*k+1].direction = -1;
                scores[j][2*k+1].score = 1.0 - deltas[didx];
            }
            std::sort(scores[j].begin(), scores[j].end());
        }
    }

    int K, L;
    std::vector<SigType> qcode;
    std::vector<Scalar> deltas;
    std::vector<std::vector<Score> > scores;
    const E2& e2lsh;
    //f :: uint32_t -> bool

    template<typename F>
    void forProbes(int idx, const F& f){
        std::priority_queue<Probe> que;
        que.emplace(0, scores[idx][0].score);

        while(!que.empty()){
            Probe cur_probe = que.top();
            que.pop();
            uint32_t code = cur_probe.getHash(idx*K, (idx+1)*K, qcode, scores[idx]);

            if(f(code)){
                return ;
            }

            int maxidx = cur_probe.maxidx;

            if(maxidx < scores[idx].size()){
                // printf("push!!\n");
                Probe shift_probe = cur_probe;
                Probe expand_probe = cur_probe;
                shift_probe.shift(scores[idx][maxidx+1].score);
                expand_probe.expand(scores[idx][maxidx+1].score);
                
                que.push(shift_probe);
                que.push(expand_probe);
            }
        }
    }
};


class MPLSH
{
public:
    //n: #points; L: #hasher; T: #probes; W: bucket width
    MPLSH(int n, int dim, int K, int L, float W) 
        : nPnts(n), dim(dim), K(K), L(L), W(W), cm(n), e2lsh(dim, K*L, W)
    {
    };
    int nPnts;
    int dim;
    int K;
    int L;
    float W;
    CountMarker cm;

    const float** data;
    E2 e2lsh;

    NDArray<2, SigType> codes;
    
    std::vector<std::unordered_map<uint32_t, std::vector<int>>> buckets;

    // void build(const std::vector<std::vector<SigType> > &codes) {
    void build(const float **data_) 
    {
        data = data_;

        // codes.reserve(n);
        printf("computing hash-code\n");
        codes.resize({ size_t(nPnts), size_t(K*L) });
        for (int i = 0; i < nPnts; ++i) {
            // codes.emplace_back(hasher->getSig(data[i]));
            e2lsh.getSig(data[i], &codes[i][0]);
        }
        printf("building index\n");

        assert(codes.lens[1] == K * L);
        SigType** codesp = codes.to_ptr();

        buckets.resize(L);
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

    void query(unsigned nCandidates, const float* query, MinK_List *list)
    {
        FancyProbeGen fpg(K, L, query, e2lsh);

        cm.clear();
        for (int j = 0; j < L; j++) {
            fpg.forProbes(j, [&](uint32_t hashcode32) -> bool{
                // printf("code=%u\n", hashcode32);
                auto it = buckets[j].find(hashcode32);
                if (it != buckets[j].end()) {
                    // printf("nCandidates=%d, left_candidates=%d\n", it->second.size(), nCandidates);
                    for (int idx : it->second) {
                        if (!cm.isMarked(idx)) {
                            float angle = calc_l2_dist(dim, data[idx], query);
                            list->insert(angle, idx + 1);
                            cm.mark(idx);
                            if (--nCandidates<=0) {
                                return true;
                            }
                        }
                    }
                }
                return false;
            });
        }
    }

    int64_t get_memory_usage()
    {
        //not implemented yet
        return 0;
    }

};