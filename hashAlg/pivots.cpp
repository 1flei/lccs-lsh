#include "pivots.h"

PivotsBinary::PivotsBinary(int d, int K, int M, int n, const Scalar** data)
    : dim(d)
    , K(K)
    , M(M)
    , nPnts(n)
    , datap(data)
{
    assert(d > 0 && K > 0 && M > 0 && n > 0 && data);
    //let nPivots*(nPivots-1)/2 > K*M;
    nPivots = (int)(sqrt(2 * K * M) + 1.5);
    pairs.reserve(nPivots * (nPivots - 1) / 2);
    pivots.resize(nPivots);
    
    for(int i=0;i<nPivots;i++){
        pivots[i] = rand()%n;
    }
    for (int i = 0; i < nPivots; i++) {
        for (int j = i + 1; j < nPivots; j++) {
            pairs.emplace_back(i, j);
        }
    }
    random_shuffle(pairs.begin(), pairs.end());
}

PivotsBinary::~PivotsBinary()
{
}

std::vector<SigType> PivotsBinary::getSig(const Scalar* data)
{
    std::vector<Scalar> dists(nPivots);
    for (int i = 0; i < nPivots; ++i) {
        dists[i] = calc_cosangle(dim, data, datap[pivots[i]]);
    }

    std::vector<SigType> ret(K);
    for (int i = 0; i < K; i++) {
        int sigi = 0;
        for (int j = 0; j < M; j++) {
            int idx = i * M + j;
            int u = pairs[idx].first;
            int v = pairs[idx].second;
            sigi = (sigi << 1) ^ (dists[u] < dists[v]);
        }
        ret[i] = sigi;
    }
    return ret;
}