#pragma once

#include <vector>
#include "../util.h"
#include <cassert>
#include <random>
#include "../Eigen/Eigen"
#include <queue>

//multi-probe version of E2Eigen
class E2MP
{
public:
    using SigType = int32_t;
    E2MP(int d, int K, double r, int maxSpan)      //dim of data object, #hasher, radius 
        :dim(d), K(K), r(r), maxSpan(maxSpan)
    {
        assert(d > 0 && K > 0);

        std::normal_distribution<double> normal(0.);
        std::uniform_real_distribution<double> uniform(0., r);
        std::random_device rd;
        std::default_random_engine rng(rd());

        p.resize(K, d);
        b.resize(K);
        for (int i = 0; i < K; i++) {
            for(int j=0;j<d;j++){
                p(i, j) = normal(rng);
            }
        }
        for (int i = 0; i < K; i++) {
            b(i) = uniform(rng) + 0.5*r;
            // b(i) = uniform(rng);
        }
        sigdim = K;

        genPertubations();
    }
    ~E2MP() {}
    std::vector<SigType> getSig(const Scalar *data)
    {
        std::vector<SigType> ret(sigdim);
        getSig(data, &ret[0]);
        return ret;
    }
    void getSig(const Scalar *data, SigType* ret)
    {
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);
        Eigen::Map<Eigen::Matrix<SigType, Eigen::Dynamic, 1> > ret_v(ret, sigdim);

        auto ret_r = (p*data_v+b)/r;
        // Eigen::Vector<Scalar, Eigen::Dynamic> ret_r = (p*data_v+b)/r + 0.5;
        ret_v = ret_r.cast<SigType>();
    }

    // static probing is better than dynamic probing
    // f :: vect<SigT>> -> last_pertubation_idx -> IO
    // template<class F> 
    // void forSig(int nProbes, const Scalar *data, const F& f)
    // {
    //     std::vector<SigType> ret(sigdim);

    //     Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);
    //     Eigen::Map<Eigen::Matrix<SigType, Eigen::Dynamic, 1> > cur_sig(&ret[0], sigdim);

    //     // Eigen::Vector<Scalar, Eigen::Dynamic> sig = (p*data_v+b)/r;
    //     auto sig = (p*data_v+b)/r;
    //     cur_sig = sig.cast<int32_t>();

    //     int probeCnt = 0;

    //     //return isEnd
    //     const auto tryProbe = [&](int last_idx) -> bool {
    //         probeCnt++;
    //         f(ret, last_idx);
    //         if(probeCnt>=nProbes){
    //             return true; 
    //         }
    //         return false;
    //     };

    //     //return isEnd
    //     const auto multi_loop = [&](int level) -> bool{
    //         auto multi_loop_impl = [&](int level, int start, auto& multi_loop_ref) -> bool {
    //             if(level==0){
    //                 return tryProbe(start);
    //             }

    //             for(int i=start;i<sigdim;i++){
    //                 ret[i]+=1;
    //                 if(multi_loop_ref(level-1, i+1, multi_loop_ref)){
    //                     return true;
    //                 }
    //                 ret[i]-=2;
    //                 if(multi_loop_ref(level-1, i+1, multi_loop_ref)){
    //                     return true;
    //                 }
    //                 ret[i]+=1;                
    //             }
    //             return false;
    //         };
    //         return multi_loop_impl(level, 0, multi_loop_impl);
    //     };

    //     //lvl 0
    //     tryProbe(-1);

    //     //lvl 1
    //     for(int i=0;i<sigdim;i++){
    //         ret[i]+=1;
    //         if(tryProbe(i)){
    //             return ;
    //         }
    //         ret[i]-=2;
    //         if(tryProbe(i)){
    //             return ;
    //         }
    //         ret[i]+=1;
    //     }

    //     //lvl 2, expand the lambda
    //     for(int i=0;i<sigdim;i++){
    //         ret[i]+=1;
    //         for(int j=i+1;j<sigdim;j++){
    //             ret[j]+=1;
    //             if(tryProbe(j)){
    //                 return ;
    //             }
    //             ret[j]-=2;
    //             if(tryProbe(j)){
    //                 return ;
    //             }
    //             ret[j]+=1;
    //         }
    //         ret[i]-=2;
    //         for(int j=i+1;j<sigdim;j++){
    //             ret[j]+=1;
    //             if(tryProbe(j)){
    //                 return ;
    //             }
    //             ret[j]-=2;
    //             if(tryProbe(j)){
    //                 return ;
    //             }
    //             ret[j]+=1;
    //         }
    //         ret[i]+=1;
    //     }
    //     //theoretically will be until 9, and then pertabute a dimension for more than one step
    //     //this may be enough for nProbe < ~O(m^9)
    //     for(int lvl=3;lvl<=9;lvl++){
    //         if(multi_loop(lvl)){
    //             return; 
    //         }
    //     }
    // }

    void genPertubations()
    {
        const int maxProbePerDim = 16;
        int maxProbe = maxProbePerDim*sigdim;

        pertubations.reserve(maxProbe+16);
        //return isend
        const auto multi_loop = [&](int level, int start, int end) -> bool{
            auto multi_loop_impl = [&](int level, int start, int end, std::vector<Pert>& cur, auto& multi_loop_ref) -> bool{
                if(level==1){
                    //will copy
                    pertubations.push_back(cur);
                    pertubations.back().emplace_back((end-1)%dim, 1);
                    pertubations.push_back(cur);
                    pertubations.back().emplace_back((end-1)%dim, -1);
                    return pertubations.size() >= maxProbe;
                }

                for(int i=start;i+level<=end;i++){
                    cur.emplace_back(i%dim, 1);
                    if(multi_loop_ref(level-1, i+1, end, cur, multi_loop_ref) ){
                        return true;
                    }
                    cur.pop_back();
                    cur.emplace_back(i%dim, -1);
                    if( multi_loop_ref(level-1, i+1, end, cur, multi_loop_ref)) {
                        return true;
                    }
                    cur.pop_back();
                }
                return false;
            };
            // pertubations.emplace_back();
            std::vector<Pert> tmp;
            tmp.emplace_back(start, 1);
            if(multi_loop_impl(level-1, start+1, end, tmp, multi_loop_impl) ){
                return true;
            }
            tmp.pop_back();
            tmp.emplace_back(start, -1);
            if(multi_loop_impl(level-1, start+1, end, tmp, multi_loop_impl) ){
                return true;
            }
            tmp.pop_back();
            return false;
        };
        //lvl 0 
        // pertubations.emplace_back();

        //lvl 1
        for(int i=0;i<sigdim;i++){
            pertubations.emplace_back();
            pertubations.back().emplace_back(i, 1);
            pertubations.emplace_back();
            pertubations.back().emplace_back(i, -1);
        }
        //lvl 2, bounded by 2**maxSpan
        [&](){
            for(int span=1;span<maxSpan;span++){
                for(int i=0;i<sigdim;i++){
                    for(int lvl=2;lvl<=span;lvl++){
                        if(multi_loop(lvl, i, i+span) ){
                            return ;
                        }
                    }
                }
            }
        }();

        printf("pertubations.size()=%d, per-dim=%d\n", pertubations.size(), pertubations.size()/sigdim);
    }

    //static probing is better than dynamic probing
    //this version considers the span of pertubation as well
    // f :: vect<SigT>> -> last_pertubation_idx -> IO
    template<class F> 
    void forSig(int nProbes, const Scalar *data, const F& f)
    {
        std::vector<SigType> ret(sigdim);

        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);
        Eigen::Map<Eigen::Matrix<SigType, Eigen::Dynamic, 1> > cur_sig(&ret[0], sigdim);

        // Eigen::Vector<Scalar, Eigen::Dynamic> sig = (p*data_v+b)/r;
        auto sig = (p*data_v+b)/r;
        cur_sig = sig.cast<int32_t>();

        f(ret, -1);

        for(int i=0;i<nProbes-1;i++){
            for(Pert& p:pertubations[i]){
                ret[p.idx] += p.shift;
            }
            f(ret, pertubations[i].back().idx);
            for(Pert& p:pertubations[i]){
                ret[p.idx] -= p.shift;
            }
        }
    }

    int64_t get_memory_usage()
    {
        return int64_t(sizeof(*this)) + (K*dim)*sizeof(Scalar) + int64_t(K)*sizeof(Scalar);
    }

    int dim, K;
    Scalar r;
    int sigdim;
    int maxSpan;
protected:
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> p;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b;

    struct Pert {
        int idx;
        int shift;
        Pert(int i, int shift): idx(i), shift(shift){}
    };

    std::vector<std::vector<Pert> > pertubations;
    
    // std::vector<Scalar> p;
    // std::vector<Scalar> b;
};

//multi-probe version of CP hasher
class CrossPolytopeMP
{
public:
    using SigType = int32_t;
};