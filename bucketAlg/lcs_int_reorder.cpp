#include "lcs_int_reorder.h"

using namespace mylccs;

LCCS_LSH_REORDER::LCCS_LSH_REORDER(int M, int step)
    :M(M), step(step)
{
    sortedLen = M*(M-1)/4 / step;
    reorderLen = sortedLen * step;
}


void LCCS_LSH_REORDER::build(NDArray<2, int32_t> &hashCodes)
{
    codep = hashCodes.to_ptr();

    // std::random_uniform_int
    std::uniform_int_distribution<int> uniform(0, M-1);
    std::random_device rd;
    std::default_random_engine rng(rd());

    for(int i=0;i<reorderLen;i++){
        reorders.push_back(uniform(rng));
        if(i > 0 && reorders[i] == reorders[i-1]){
            reorders[i] = (reorders[i] + 1)%M;
        }
    }

    nPnts = hashCodes.lens[0];

    checked.resize(nPnts);

    init_sorted_idx();
    init_next_link();
}


std::tuple<int, bool> LCCS_LSH_REORDER::match_util(int pidxx, int pidxy, int dimIdxBeg, int match)
{
    //to make sure the order is correct
    if(match>0){
        --match;
    }
    int matchloc = dimIdxBeg+match;
    int restlen = reorderLen - match;
    for (int i = 0; i < restlen; i++) {
        int matchDimIdx = (matchloc + i);
        int x_i = get_hash_code(pidxx, matchDimIdx);
        int y_i = get_hash_code(pidxy, matchDimIdx);
        if (x_i!=y_i) return std::make_tuple(match, x_i < y_i);
        match++;
    }
    //equal
    match = reorderLen;
    return std::make_tuple(match, false);
} 


std::tuple<int, bool> LCCS_LSH_REORDER::match_util(const int32_t* x, const int32_t* y, int dimIdxBeg, int match)
{
    //to make sure the order is correct
    if(match>0){
        --match;
    }
    int matchloc = dimIdxBeg+match;
    int restlen = reorderLen - match;
    for (int i = 0; i < restlen; i++) {
        int matchDimIdx = (matchloc + i);
        int x_i = get_hash_code(x, matchDimIdx);
        int y_i = get_hash_code(y, matchDimIdx);
        if (x_i!=y_i) return std::make_tuple(match, x_i < y_i);
        match++;
    }
    //equal
    match = reorderLen;
    return std::make_tuple(match, false);
} 





//lowidx, lowlen, highlen
//assume the length of matching of low and high are lowlen and
std::tuple<int, int, int> LCCS_LSH_REORDER::get_loc_scan(const int32_t* q, int searchLoc, int low, int lowlen, int high, int highlen)
{
    int lastlen = lowlen;
    int minlen = std::min(lowlen, highlen);
    for(int i=low+1;i<high;i++){
        const int32_t* dpi = codep[sorted_idx[searchLoc][i]];
        auto [ilen, iless] = match_util(q, dpi, searchLoc*step, minlen);

        if(iless){
            return std::make_tuple(i-1, lastlen, ilen);
        }
        lastlen = ilen;
    }
    //reach the end
    return std::make_tuple(high-1, lastlen, highlen);
}

//binary search with linear when interval is small
//return  (lowidx, lowlen, highlen) that q \in [lowidx, highidx) or lowidx==0 or highidx==N-1
std::tuple<int, int, int> LCCS_LSH_REORDER::get_loc_mixed(const int32_t* query, int searchloc, int low, int lowlen, int high, int highlen)
{
    static const int SCAN_SIZE = 4;
    // const auto& sorted_idx_qloc = sorted_idx[searchloc];

    //return datap->at(idx[searchloc][a]) < query
    const auto& cmp = [&](int a, int match=0) -> std::tuple<int, bool> {
        const int32_t* dpa = codep[sorted_idx[searchloc][a]];
        return match_util(query, dpa, searchloc*step, match);
    };


    while (low < high - SCAN_SIZE) {
        int curlen = std::min(lowlen, highlen);
        //binary search
        int mid = (low + high) / 2;
        auto [midlen, midisless] = cmp(mid, curlen);
        // f(searchloc, mid, midlen, midisless);
        
        //if q<mid:
        if(midisless) {
            //which means mid < query
            high = mid;
            highlen = midlen;
        } else{
            low = mid;
            lowlen = midlen;
        }
    }

    return get_loc_scan(query, searchloc, low, lowlen, high, highlen);
}

//return  (lowidx, lowlen, highlen) that q \in [lowidx, highidx) or lowidx==0 or highidx==N-1
std::tuple<int, int, int> LCCS_LSH_REORDER::get_loc(const int32_t* query, int searchloc)
{
    int low = 0;
    int high = nPnts-1;
    const int32_t* dplow = codep[sorted_idx[searchloc][low]];
    const int32_t* dphigh = codep[sorted_idx[searchloc][high]];
    auto [lowlen, lowisless] = match_util(query, dplow, searchloc*step, 0);
    auto [highlen, highisless] = match_util(query, dphigh, searchloc*step, 0);

    return get_loc_mixed(query, searchloc, low, lowlen, high, highlen);

}




void LCCS_LSH_REORDER::init_sorted_idx()
{
    sorted_idx.resize(sortedLen);
    for (int d = 0; d < sortedLen; d++) {
        sorted_idx[d].resize(nPnts);
        for (int i = 0; i < nPnts; i++) {
            sorted_idx[d][i] = i;
        }

        int matchDimIdx = d*step;
        
        const auto& cmpf = [matchDimIdx, this](int32_t a, int32_t b) -> bool {
            auto [len, isLess] = match_util(a, b, matchDimIdx, 0); 
            return isLess;
        };

        sort(sorted_idx[d].begin(), sorted_idx[d].end(), cmpf);
    }
}

void LCCS_LSH_REORDER::init_next_link()
{
    std::unordered_map<int, int> nextMap;
    
    next_link.resize(sortedLen);
    for(int d=sortedLen-1;d>=0;--d){
        nextMap.clear();

        next_link[d].resize(nPnts);
        int nextd = (d+1)%sortedLen;
        for(int i=0;i<nPnts;i++){
            nextMap[sorted_idx[nextd][i]] = i;
        }

        for(int i=0;i<nPnts;i++){
            int dataidxi = sorted_idx[d][i];
            next_link[d][i] = nextMap[dataidxi];
        }
    }
}