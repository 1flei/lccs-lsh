#pragma once

#include "../hashAlg/e2eigen.h"
#include "../hashAlg/polytope.h"
#include <queue>

//multi-probe version of E2Eigen
class E2MP
{
public:
    using SigType = int32_t;
    E2MP(int d, int K, double r, int maxProbePerDim=16)      //dim of data object, #hasher, radius 
        :dim(d), K(K), r(r), maxProbePerDim(maxProbePerDim)
    {
        assert(d > 0 && K > 0);

        std::normal_distribution<double> normal(0.);
        std::uniform_real_distribution<double> uniform(0., r);
        // std::random_device rd;
        // std::default_random_engine rng(rd());
        std::default_random_engine rng(GLOBAL_SEED);

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
        probe_indices.resize(probes.size());
        for(int i=0;i<probe_indices.size();i++){
            probe_indices[i] = i;
        }
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

    void genPertubations()
    {
        int maxProbe = 128*sigdim;

        probes.reserve(maxProbe+16);
        //return isend
        //__P__
        for(int i=0;i<sigdim;i++){
            probes.emplace_back(9, i, 1);
        }
        for(int i=0;i<sigdim;i++){
            probes.emplace_back(9, i, -1);
        }
        if(probes.size() > maxProbe){
            return ;
        }

        //__PP__
        for(int i=0;i<sigdim;i++){
            probes.emplace_back(18, i, 1,  (i+1)%sigdim, 1);
            probes.emplace_back(18, i, 1,  (i+1)%sigdim, -1);
            probes.emplace_back(18, i, -1, (i+1)%sigdim, 1);
            probes.emplace_back(18, i, -1, (i+1)%sigdim, -1);
        }
        if(probes.size() > maxProbe){
            return ;
        }

        //__P_P__
        for(int i=0;i<sigdim;i++){
            probes.emplace_back(23, i, 1,  (i+2)%sigdim, 1);
            probes.emplace_back(23, i, 1,  (i+2)%sigdim, -1);
            probes.emplace_back(23, i, -1, (i+2)%sigdim, 1);
            probes.emplace_back(23, i, -1, (i+2)%sigdim, -1);
        }
        if(probes.size() > maxProbe){
            return ;
        }

        //__P__P__
        for(int i=0;i<sigdim;i++){
            probes.emplace_back(28, i, 1,  (i+3)%sigdim, 1);
            probes.emplace_back(28, i, 1,  (i+3)%sigdim, -1);
            probes.emplace_back(28, i, -1, (i+3)%sigdim, 1);
            probes.emplace_back(28, i, -1, (i+3)%sigdim, -1);
        }
        if(probes.size() > maxProbe){
            return ;
        }

        //__PPP__
        for(int i=0;i<sigdim;i++){
            probes.emplace_back(27, i, 1,  (i+1)%sigdim, 1,  (i+2)%sigdim, 1);
            probes.emplace_back(27, i, 1,  (i+1)%sigdim, 1,  (i+2)%sigdim, -1);
            probes.emplace_back(27, i, 1,  (i+1)%sigdim, -1, (i+2)%sigdim, 1);
            probes.emplace_back(27, i, 1,  (i+1)%sigdim, -1, (i+2)%sigdim, -1);
            probes.emplace_back(27, i, -1, (i+1)%sigdim, 1,  (i+2)%sigdim, 1);
            probes.emplace_back(27, i, -1, (i+1)%sigdim, 1,  (i+2)%sigdim, -1);
            probes.emplace_back(27, i, -1, (i+1)%sigdim, -1, (i+2)%sigdim, 1);
            probes.emplace_back(27, i, -1, (i+1)%sigdim, -1, (i+2)%sigdim, -1);
        }
        if(probes.size() > maxProbe){
            return ;
        }
    }

    //static probing 
    //this version considers the span of pertubation as well
    // f :: vect<SigT>> -> last_pertubation_idx -> IO
    template<class F> 
    void forSigStatic(int nProbes, const Scalar *data, const F& f)
    {
        std::vector<SigType> ret(sigdim);

        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);
        Eigen::Map<Eigen::Matrix<SigType, Eigen::Dynamic, 1> > cur_sig(&ret[0], sigdim);

        // Eigen::Vector<Scalar, Eigen::Dynamic> sig = (p*data_v+b)/r;
        auto sig = (p*data_v+b)/r;
        cur_sig = sig.cast<int32_t>();

        f(ret, -1);

        for(int i=0;i<nProbes-1;i++){
            probes[i].for_pertubated(ret, f);
        }
    }

    struct Probe
    {
        static const int MAX_PERTS = 3;
        double score;
        std::array<int32_t, MAX_PERTS> perts;
        std::array<int8_t, MAX_PERTS> dirs;

        template<class FCode>
        void for_pertubated(std::vector<SigType>& codes, const FCode& f){
            // printf("score= %f\n", score);
            // printf("pert=  %d, %d, %d\n", perts[0], perts[1], perts[2]);
            // printf("dirs=  %d, %d, %d\n", dirs[0], dirs[1], dirs[2]);
            for(int i=0;i<MAX_PERTS;i++){
                int p = perts[i];
                int dir = dirs[i];
                codes[p] += dir;
            }
            f(codes, perts[0]);
            for(int i=0;i<MAX_PERTS;i++){
                int p = perts[i];
                int dir = dirs[i];
                codes[p] -= dir;
            }
        }

        bool operator<(const Probe& p) const{
            return score < p.score;
        }

        void set_score(double score_) {
            score = score_;
        }

        Probe(double score, int p0, int8_t shift0)
            : score(score), perts{{p0}}, dirs{{shift0, 0, 0}}
        {} 
        Probe(double score, int p0, int8_t shift0, int p1, int8_t shift1)
            : score(score), perts{{p0, p1}}, dirs{{shift0, shift1, 0}}
        {} 
        Probe(double score, int p0, int8_t shift0, int p1, int8_t shift1, int p2, int8_t shift2)
            : score(score), perts{{p0, p1, p2}}, dirs{{shift0, shift1, shift2}}
        {} 
    };

    //dynamic probing
    template<class F> 
    void forSig(int nProbes, const Scalar *data, const F& f)
    {
        std::vector<SigType> ret(sigdim);

        using VectF = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        using VectI = Eigen::Matrix<SigType, Eigen::Dynamic, 1>;

        Eigen::Map<const VectF> data_v(data, dim);
        Eigen::Map<VectI > cur_sig(&ret[0], sigdim);

        // Eigen::Vector<Scalar, Eigen::Dynamic> sig = (p*data_v+b)/r;
        VectF sig = (p*data_v+b)/r;
        VectF sig_floor = sig.array().floor();
        cur_sig = sig.cast<int32_t>();
        
        f(ret, -1);
        if(nProbes<=1){
            return ;
        }
        --nProbes;

        // auto diffp = (cur_sig  -sig) +1.5;
        // auto diffn = (sig - cur_sig) +0.5;
        VectF diffp = (sig_floor  -sig);
        // VectF diffn = (sig - sig_floor);

        auto updateProbe = [&](Probe& p){
            double score = 0.;
            for(int i=0;i<Probe::MAX_PERTS;i++){
                if(p.dirs[i] < 0){
                    double sqrts = -diffp(p.perts[i])+0.5;
                    score += sqrts*sqrts;
                } else if(p.dirs[i] >0){
                    double sqrts = diffp(p.perts[i])+1.5;
                    score += sqrts*sqrts;
                }
            }
            p.set_score(score);
        };

        for(int i=0;i<probes.size();i++){
            updateProbe(probes[i]);
        }
        std::nth_element(probe_indices.begin(), probe_indices.begin()+nProbes, probe_indices.end(), [&](int a, int b){
            return probes[a] < probes[b];
        });
        // std::sort(probe_indices.begin(), probe_indices.begin()+nProbes, [&](int a, int b){
        //     return probes[a] < probes[b];
        // });
        for(int i=0;i<nProbes-1;i++){
            // printf("probe-%d\n", i);
            probes[probe_indices[i]].for_pertubated(ret, f);
        }
    }


    int64_t get_memory_usage()
    {
        return int64_t(sizeof(*this)) + (K*dim)*sizeof(Scalar) + int64_t(K)*sizeof(Scalar);
    }

    int dim, K;
    Scalar r;
    int sigdim;
    int maxProbePerDim;
protected:
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> p;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b;

    struct Pert {
        int idx;
        int shift;
        Pert(int i, int shift): idx(i), shift(shift){}
    };

    // std::vector<std::vector<Pert> > pertubations;
    std::vector<Probe> probes;
    std::vector<int> probe_indices;
    
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
    int sigdim;
    int num_rotations;
    int last_cp_dim;
    int maxProbePerDim;
    Hasher hasher;
    Transformer transfromer;
    TransformedVectorType transformedVec;

    CrossPolytopeMP(int vector_dim,
                int l, int maxProbePerDim=8, int num_rotations=1,
                int seed=GLOBAL_SEED)
        :dim(vector_dim), 
         sigdim(l), 
         num_rotations(num_rotations), 
         last_cp_dim(vector_dim), 
         //k=1
         maxProbePerDim(maxProbePerDim), 
         hasher(vector_dim, 1, l, num_rotations, last_cp_dim, seed), 
         transfromer(hasher)
    {
        hasher.reserve_transformed_vector_memory(&transformedVec);
        printf("cp_mp_hasher: dim=%d, l=%d, maxProbePerDim=%d\n", vector_dim, l, maxProbePerDim);

        idx.resize(sigdim);
        for(int i=0;i<sigdim;i++){
            idx[i].resize(vector_dim*2);
            for(int j=0;j<vector_dim*2;j++){
                idx[i][j] = j;
            }
        }
    }

    std::vector<SigType> getSig(const Scalar *data)
    {
        std::vector<SigType> ret(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);

        hasher.hash(data_v, &ret, &transformedVec);
        return ret;
    }

    void getSig(const Scalar *data, SigType* ret)
    {
        std::vector<SigType> res(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);

        hasher.hash(data_v, &res, &transformedVec);
        std::copy(res.begin(), res.end(), ret);
    }

    struct ProbeCandidate
    {
        double hashValDiff;
        int loc;
        int value;
        ProbeCandidate(double d, int loc, int v): hashValDiff(d), loc(loc), value(v) {}

        bool operator<(const ProbeCandidate& p) const {
            return hashValDiff < p.hashValDiff;
        }
    };

    //cannot easily apply static pertubations, will do dynamic one
    //this version considers the span of pertubation as well
    // f :: vect<SigT>> -> last_pertubation_idx -> IO
    template<class F> 
    void forSig(int nProbes, const Scalar *data, const F& f)
    {
        if(nProbes <= 0){
            return ;
        }

        std::vector<SigType> ret(sigdim);
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > data_v(data, dim);
        if(nProbes <= 1){
            hasher.hash(data_v, &ret, &transformedVec);
            f(ret, -1);
            return ;
        }

        std::unique_ptr<HashHelper> fht = std::make_unique<HashHelper>(dim);
        hasher.compute_rotated_vectors(data_v, &transformedVec, fht.get());
        // const int maxProbePerDim = 16;
        
        for(int i=0;i<sigdim;i++){
            const auto cmpf = [&](int a, int b){
                double x = a<dim ? transformedVec[i][a] : -transformedVec[i][a-dim];
                double y = b<dim ? transformedVec[i][b] : -transformedVec[i][b-dim];
                return x > y;
            };
            std::nth_element(idx[i].begin(), idx[i].begin()+maxProbePerDim+1, idx[i].end(), cmpf);
            std::sort(idx[i].begin(), idx[i].begin()+maxProbePerDim+1, cmpf);

            // printVec(&transformedVec[i][0], dim);
            // printf("idx=%d, %d, ... \n\n", idx[i][0], idx[i][1]);

            ret[i] = idx[i][0];
        }

        //first probe
        f(ret, -1);
        int probeCnt = 1;
        for(int j=1;j<maxProbePerDim+1;j++){
            for(int i=0;i<idx.size();i++){
                int tmp = ret[i];
                ret[i] = idx[i][j];
                f(ret, i);
                ret[i] = tmp;
            }
        }
    }

    std::vector<std::vector<int> > idx;

    int64_t get_memory_usage()
    {
        return 0;
    }
};