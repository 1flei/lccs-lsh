#pragma once

#include "../hashAlg/e2eigen.h"
#include "../hashAlg/polytope.h"

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
            f(ret, pertubations[i][0].idx);
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

    // using SigType = uint32_t;
    // using Hasher = falconn::core::CrossPolytopeHashBase<falconn::core::CrossPolytopeHashDense<Scalar, SigType>, std::vector<float>, Scalar, SigType>;
    using Hasher = falconn::core::CrossPolytopeHashDense<Scalar, SigType>;
    using Transformer = Hasher::HashTransformation;
    using TransformedVectorType = Hasher::TransformedVectorType;
    using HashHelper = falconn::core::cp_hash_helpers::FHTHelper<Scalar>;

    int dim;
    int l;
    int num_rotations;
    int last_cp_dim;
    int maxSpan;
    int sigdim;
    Hasher hasher;
    Transformer transfromer;
    TransformedVectorType transformedVec;

    CrossPolytopeMP(int vector_dim,
                int l, int last_cp_dim=16, int maxSpan=8, int num_rotations=1,
                int seed=666)
        :dim(vector_dim), 
         sigdim(l), 
         num_rotations(num_rotations), 
         last_cp_dim(last_cp_dim), 
         //k=1
         maxSpan(maxSpan), 
         hasher(vector_dim, 1, l, num_rotations, last_cp_dim, seed), 
         transfromer(hasher)
    {
        hasher.reserve_transformed_vector_memory(&transformedVec);
    }

    std::vector<SigType> getSig(const Scalar *data)
    {
        std::vector<SigType> ret(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);

        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> normalized_data_v = data_v.normalized();
        hasher.hash(normalized_data_v, &ret, &transformedVec);
        return ret;
    }

    void getSig(const Scalar *data, SigType* ret)
    {
        std::vector<SigType> res(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);

        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> normalized_data_v = data_v.normalized();
        hasher.hash(normalized_data_v, &res, &transformedVec);
        std::copy(res.begin(), res.end(), ret);
    }


    //cannot easily apply static pertubations, will do dynamic one
    //this version considers the span of pertubation as well
    // f :: vect<SigT>> -> last_pertubation_idx -> IO
    template<class F> 
    void forSig(int nProbes, const Scalar *data, const F& f)
    {
        std::vector<SigType> ret(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);

        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> normalized_data_v = data_v.normalized();
        std::unique_ptr<HashHelper> fht = std::make_unique<HashHelper>(dim);
        hasher.compute_rotated_vectors(normalized_data_v, &transformedVec, fht.get());
        // hasher.compute_cp_hashes(transformedVec, 1, l_, rotation_dim_,
        //               log_rotation_dim_, last_cp_dim_, last_cp_log_dim_,
        //               result);
        printf("transformed vector\n");
        for(int i=0;i<transformedVec.size();i++){
            printVec(&transformedVec[i][0], dim);
        }

        hasher.hash(normalized_data_v, &ret, &transformedVec);

        f(ret, -1);

        // for(int i=0;i<nProbes-1;i++){
        //     for(Pert& p:pertubations[i]){
        //         ret[p.idx] += p.shift;
        //     }
        //     f(ret, pertubations[i][0].idx);
        //     for(Pert& p:pertubations[i]){
        //         ret[p.idx] -= p.shift;
        //     }
        // }
    }

    int64_t get_memory_usage()
    {
        return 0;
    }
};