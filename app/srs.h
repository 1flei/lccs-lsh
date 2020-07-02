#pragma once


#include "../SRS/SRSInMemory.h"
#include "../util.h"
#include "../myndarray.h"
#include <unordered_map>
#include "../register.h"
#include "../benchmark_util.h"
#include <queue>
#include "memory.h"
#include <boost/math/distributions/chi_squared.hpp>
#include "../Eigen/Eigen"


struct MYSRS
{
    MYSRS(int nPnts, int dim, int proj_dim) 
        :nPnts(nPnts), dim(dim), proj_dim(proj_dim)
    {
        std::normal_distribution<double> normal(0.);
        // std::random_device rd;
        // std::default_random_engine rng(rd());
        std::default_random_engine rng(GLOBAL_SEED);

        p.resize(proj_dim, dim);
        for (int i = 0; i < proj_dim; i++) {
            for(int j=0;j<dim;j++){
                p(i, j) = normal(rng);
            }
        }
    }

    void getSig(const float *data, float* ret)
    {
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, 1> > data_v(data, dim);
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1> > ret_v(ret, proj_dim);

        ret_v = p*data_v;
    }

    // void build(const std::vector<std::vector<SigType> > &codes) {
    void build(const float **data_) 
    {
        data = (float**)data_;
        proj_data.reset(new Proj_data(nPnts, proj_dim, nullptr));
        //SRS_cover_tree requires continuous data
        for (int i = 0; i < nPnts; ++i)
        {
            //use 
            float* proj_datai = &proj_data->data[i * proj_dim];
            getSig(data_[i], proj_datai);
        }
        //now proj_data will have all projected data
        tree_index.reset(new SRS_Cover_Tree(nPnts, proj_dim, proj_data.get()) );
    }

    void query(unsigned nCandidates, const float* query, MinK_List *list)
    {
        std::vector<float> q_proj_tmp(proj_dim);
        // float *q_proj = new float[proj_dim];
        float* q_proj = &q_proj_tmp[0];
        getSig(query, q_proj);

        tree_index->init_search(q_proj);
        for(int count=0; count < nCandidates; count++) {
            res_pair cover_tree_res = tree_index->increm_knn_search_compressed();
            double dist = calc_l2_dist(dim, data[cover_tree_res.id], query);

            // double dist_proj = calc_l2_dist(proj_dim, &proj_data->data[cover_tree_res.id*proj_dim], q_proj);
            // printf("dist=%f, res.dist=%f\n", dist, cover_tree_res.dist);

            
            list->insert(dist, cover_tree_res.id+1);
        }

        tree_index->finish_search();
        return;
    }

    std::unique_ptr<Proj_data> proj_data;
    std::unique_ptr<SRS_Cover_Tree> tree_index;
    float **data;

    int nPnts, dim, proj_dim;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> p;
};