#pragma once

#include <vector>
#include "../def.h"
#include "../util.h"

//
class PivotsBinary
{
public:
    PivotsBinary(int d, int K, int M, int n, const Scalar **data);      //dim of data object, #hasher, #bits per hasher, pointer to dataset 
    ~PivotsBinary();
    std::vector<SigType> getSig(const Scalar *data);
protected:
    int dim, K, M, nPnts;

    int nPivots;
    std::vector<int> pivots;
    std::vector<std::pair<int, int> > pairs;
    const Scalar** datap;
};