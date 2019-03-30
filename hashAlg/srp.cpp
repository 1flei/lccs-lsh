#include "srp.h"

using namespace std;
using namespace MyCallbackRegister;

bool SRP::registered = registerCallback("srp", "msg", [&](){
    cout << "srp registered" << endl;
    cout << algAs<string>("msg") << endl;
});

SRP::SRP(int d, int K):dim(d), M(K), p(K*d) 
{
    assert(d>0 && K>0);

    std::normal_distribution<double> normalDistribution(0.);
    std::random_device rd;
    std::default_random_engine rng(rd() );
    for(int i=0;i<K*d;i++){
        p[i] = normalDistribution(rng);
    }
}

SRP::~SRP()
{
}

vector<bool> SRP::getSig(Scalar* data)
{
    vector<bool> ret(M);
    for(int i=0;i<M;i++){
        ret[i] = calc_inner_product(dim, &p[i*dim], data);
    }
    return ret;
}