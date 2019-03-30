#pragma once
#include <vector>
#include "../util.h"
#include <cassert>
#include <random>
#include "../register.h"

//Random Sign Projection
class SRP
{
public:
    SRP(int d, int K);      //dim of data object, #hasher 
    ~SRP();

    //
    template<typename Container>
    std::vector<bool> getSig(Container& data){
        std::vector<bool> ret(M);

        for(int i=0;i<M;i++){
            Scalar si = 0.;
            int j = 0;
            for(Scalar dj:data){
                int idx = i*dim+j;
                si += p[idx]*dj;

                if(++j>=dim){
                    break;
                }
            }
            ret[i] = si;
        }
    }
    
    template<typename Containers>
    std::vector<std::vector<bool> > getSigs(Containers& data){
        std::vector<std::vector<bool> > ret;
        for(auto& datai:data){
            //push rvalue
            ret.emplace_back(getSig(datai));
        }
        return ret;
    }

    //unsafe!
    std::vector<bool> getSig(Scalar* data);

    std::vector<std::vector<bool> > getSigs(Scalar** data);

protected:
    std::vector<Scalar> p;
    int dim, M;

    static bool registered;
};
