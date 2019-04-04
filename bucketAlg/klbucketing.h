#pragma once

#include "../util.h"
#include <unordered_map>

//bucketing algorithm using (k,L)-setting
//i.e., verify as per one match, from this point of view, this is like L_-inf bucketing
class KLBucketing
{
public:
    KLBucketing(int d) : dim(d), buckets(d), cm(0){};
    int dim;
    const std::vector<std::vector<SigType> > *codesp;

    std::vector<std::unordered_map<SigType, std::vector<int> > > buckets;

    void build(const std::vector<std::vector<SigType> > &codes) {
        //build hash table based one codes
        for(int i=0;i<codes.size();i++){
            for(int j=0;j<dim;j++){
                buckets[j][codes[i][j]].push_back(i);
            }
        }

        cm.resize(codes.size());
    }
    
    template<typename F>
    void forCandidates(int nCandidates, const std::vector<SigType> &qcode, const F& f) {
        cm.clear();
        for(int j=0;j<dim;j++){
            if(auto it = buckets[j].find(qcode[j]); it !=buckets[j].end()){
                for(int idx:it->second){
                    if(!cm.isMarked(idx)){
                        f(idx);
                        cm.mark(idx);
                        if(!--nCandidates){
                            return ;
                        }
                    }
                }
            }
        }
    }

    CountMarker cm;
};