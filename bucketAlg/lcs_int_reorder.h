#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <queue>
#include <set>
#include <cinttypes>
#include <array>
#include <functional>
#include "../myTimer.h"
#include <unordered_map>
#include "../myndarray.h"
#include "../util.h"

namespace mylcs 
{

	class LCCS_LSH_REORDER
	{
	public:
		LCCS_LSH_REORDER(int M, int step);
        void build(NDArray<2, int32_t> &hashCodes);

		inline int32_t get_hash_code(int pntidx, int dimidx)
		{
			return codep[pntidx][reorders[dimidx%reorderLen]];
		}
		inline int32_t get_hash_code(const int32_t* rawp, int dimidx)
		{
			return rawp[reorders[dimidx%reorderLen]];
		}

		//match_util :: [sigs] a-> [sigs] b->int loc-> (int, bool)
		//given two indices pidxx, pidxy in codep
        //return matchlen and bool(x<y)
		std::tuple<int, bool> match_util(int pidxx, int pidxy, int dimIdxBeg, int match);

		//match_util :: [sigs] a-> [sigs] b->int loc-> (int, bool)
		//given two raw pointer x and y
        //return matchlen and bool(x<y)
		std::tuple<int, bool> match_util(const int32_t* x, const int32_t* y, int loc, int match);

		//lowidx, lowlen, highlen
		//assume the length of matching of low and high are lowlen and highlen
		std::tuple<int, int, int> get_loc_scan(const int32_t* query, int loc, int low, int lowlen, int high, int highlen);

		//binary search with linear search when interval is small
		//return  (lowidx, lowlen, highlen) that q \in [lowidx, highidx) or lowidx==0 or highidx==N-1
		std::tuple<int, int, int> get_loc_mixed(const int32_t* query, int qloc, int low, int lowlen, int high, int highlen);

		//return  (lowidx, lowlen, highlen) that q \in [lowidx, highidx) or lowidx==0 or highidx==N-1
		std::tuple<int, int, int> get_loc(const int32_t* query, int qloc);


		void init_sorted_idx();
		void init_next_link();


		int64_t get_memory_usage()
		{
			return 0;
		}

		template<typename F>
		void for_candidates(int scanStep, const std::vector<int32_t>& query, const F& f) 
		{
			const int32_t* queryp = (const int32_t*)&query[0];
			
            checked.clear();
            const auto& tryCheck = [&](int dataidx){
				if (!checked.isMarked(dataidx)) {
                    checked.mark(dataidx);
					f(dataidx);
				}
            };

            const auto& tryCheckLoc = [&](int curidx, int d){
				// MyTimer::pusht();
				for(int i=curidx;i>=0 && curidx-i<scanStep;--i){
					int matchingIdx = sorted_idx[d][i];
					tryCheck(matchingIdx);
				}
				for(int i=curidx+1;i<nPnts && i-curidx-1<scanStep;i++){
					int matchingIdx = sorted_idx[d][i];
					tryCheck(matchingIdx);
				}
				// checktime += MyTimer::popt();
            };

			// MyTimer::pusht();
			auto [curidx, lowlen, highlen] = get_loc(queryp, 0);
			// searchtime += MyTimer::popt();

			tryCheckLoc(curidx, 0);
			// printf("cur, lowlen, highlen: %d, %d, %d, [lowidx=%d, highidx=%d]\n", curidx, lowlen, highlen, sorted_idx[0][curidx], sorted_idx[0][curidx+1]);
			
			for(int d=1;d<sortedLen;d++){
				// MyTimer::pusht();
				int lowidx = next_link[d-1][curidx];
				int highidx = next_link[d-1][curidx+1];

				if(lowlen < step || highlen < step){
					std::tie(curidx, lowlen, highlen) = get_loc(queryp, d);
				} else{
					if(lowlen!=reorderLen){
						lowlen -= step;
					}
					if(highlen!=reorderLen){
						highlen -= step;
					}
					// std::tie(curidx, lowlen, highlen) = get_loc_scan(queryp, d, lowidx, lowlen, highidx, highlen);
					std::tie(curidx, lowlen, highlen) = get_loc_mixed(queryp, d, lowidx, lowlen, highidx, highlen);
				}
				// printf("cur, lowlen, highlen: %d, %d, %d, [lowidx=%d, highidx=%d]\n", curidx, lowlen, highlen, sorted_idx[d][curidx], sorted_idx[d][curidx+1]);
				// searchtime += MyTimer::popt();
				tryCheckLoc(curidx, d);
			}

			// printf("checktime=%f, searchtime=%f\n", checktime, searchtime);
		}


		int M;
		int step;
		int reorderLen;
		int sortedLen;
		int nPnts;
		std::vector<int> reorders;
		int32_t** codep;
        CountMarker checked;

		std::vector<std::vector<int32_t> > sorted_idx;
		std::vector<std::vector<int32_t> > next_link;

	};

}