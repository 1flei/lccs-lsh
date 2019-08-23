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

namespace mylcs 
{//LCS by sort
	class LCS_SORT_INT
	{
	public:
		LCS_SORT_INT(int dim, int step)
            :dim(dim), step(step)
		{
			// assert(dim%64==0);
			nSearchLoc = ((dim-1)/step)+1;
		}

        void build(NDArray<2, int32_t> &data)
        {
            datap = data.to_ptr();
			// nPnts = datap->size();
			nPnts = data.lens[0];

            checked.resize(nPnts);

			// printf("init!!\n");

			init_sorted_idx();
			init_next_link();
        }

		struct DataIdx
		{
			int32_t idx;

			DataIdx() {}
			DataIdx(int32_t i): idx(i){}
		};

		int dim;
		int nPnts;
		// double alpha;
		int step;
		int nSearchLoc;
		int32_t** datap;
		// std::vector<std::vector<int32_t> > sorted_idx;
		std::vector<std::vector<DataIdx> > sorted_idx;
		std::vector<std::vector<int32_t> > next_link;

		int logn;

        CountMarker checked;
        // CountMarkerU<uint16_t> inque;

		inline const int32_t* get_datap(int qloc, int idx)
		{
			return datap[sorted_idx[qloc][idx].idx];
		}
		inline const int32_t* get_datap(const DataIdx& idx)
		{
			// return &datap->at(idx.idx)[0];
			return datap[(idx.idx)];
		}
		inline int getidx(const DataIdx& idx)
		{
			return idx.idx;
		}
		inline int getidx(int qloc, int idx)
		{
			return sorted_idx[qloc][idx].idx;
		}

		inline int dim64()
		{
			return dim / 64;
		}

		//bucket-sort leveraging extents
		void init_sorted_idx()
		{
			int TableSize = (1<<logn)+1;
			sorted_idx.resize(dim);
			for (int d = 0; d < dim; d++) {
				sorted_idx[d].resize(nPnts);
				for (int i = 0; i < nPnts; i++) {
					sorted_idx[d][i] = {i};
				}
				
				const auto& cmpf = [d, this](const DataIdx& a, const DataIdx& b) -> bool {
					auto [len, isLess] = match_util(a, b, d, 0); 
					return isLess;
				};

				sort(sorted_idx[d].begin(), sorted_idx[d].end(), cmpf);
			}
		}

		void init_next_link()
		{
			std::unordered_map<int, int> nextMap;
			
			next_link.resize(dim);
			for(int d=dim-1;d>=0;--d){
				nextMap.clear();

				next_link[d].resize(nPnts);
				int nextd = (d+step)%dim;
				for(int i=0;i<nPnts;i++){
					// nextMap[sorted_idx[nextd][i]] = i;
					nextMap[getidx(nextd, i)] = i;
				}

				for(int i=0;i<nPnts;i++){
					// int dataidxi = sorted_idx[d][i];
					int dataidxi = getidx(d, i);
					next_link[d][i] = nextMap[dataidxi];
				}
			}
		}

		//lowidx, lowlen, highlen
		//assume the length of matching of low and high are lowlen and
		std::tuple<int, int, int> get_loc_scan(const int32_t* q, int loc, int low, int lowlen, int high, int highlen)
		{
			int lastlen = lowlen;
			int minlen = std::min(lowlen, highlen);
			for(int i=low+1;i<high;i++){
				const int32_t* dpi = get_datap(loc, i);
				auto [ilen, iless] = match_util(q, dpi, loc, minlen);

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
		std::tuple<int, int, int> get_loc_mixed(const int32_t* query, int qloc, int low, int lowlen, int high, int highlen)
		{
			static const int SCAN_SIZE = 4;
			// const auto& sorted_idx_qloc = sorted_idx[qloc];

			//return datap->at(idx[qloc][a]) < query
			const auto& cmp = [&](int a, int match=0) -> std::tuple<int, bool> {
				const int32_t* dpa = get_datap(qloc, a);
				return match_util(query, dpa, qloc, match);
			};


			while (low < high - SCAN_SIZE) {
				int curlen = std::min(lowlen, highlen);
				//binary search
				int mid = (low + high) / 2;
				auto [midlen, midisless] = cmp(mid, curlen);
				// f(qloc, mid, midlen, midisless);
				
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

			return get_loc_scan(query, qloc, low, lowlen, high, highlen);
		}

		//return  (lowidx, lowlen, highlen) that q \in [lowidx, highidx) or lowidx==0 or highidx==N-1
		std::tuple<int, int, int> get_loc(const int32_t* query, int qloc)
		{
			int low = 0;
			int high = nPnts-1;
			const int32_t* dplow = get_datap(qloc, low);
			const int32_t* dphigh = get_datap(qloc, high);
			auto [lowlen, lowisless] = match_util(query, dplow, qloc, 0);
			auto [highlen, highisless] = match_util(query, dphigh, qloc, 0);

			return get_loc_mixed(query, qloc, low, lowlen, high, highlen);

		}


		//match_util :: [sigs] a-> [sigs] b->int loc-> (int, bool)
		//return matchlen and bool(x<y)
		std::tuple<int, bool> match_util(const int32_t* x, const int32_t* y, int loc, int match)
		{
			//to make sure the order is correct
			if(match>0){
				--match;
			}
			int matchloc = (loc+match)%dim;
			int restlen = dim - match;
			for (int i = 0; i < restlen; i++) {
				int idx = (matchloc + i) % dim;
				int x_i = x[idx];
				int y_i = y[idx];
				if (x_i!=y_i) return std::make_tuple(match, x_i < y_i);
				match++;
			}
			//equal
			match = dim;
			return std::make_tuple(match, false);
		} 

		std::tuple<int, bool> match_util(const DataIdx& x, const DataIdx& y, int loc, int match)
		{
			//to make sure the order is correct
			const int32_t* xp = datap[x.idx];
			const int32_t* yp = datap[y.idx];
			return match_util(xp, yp, loc, match);
		} 

		struct CandidateLoc {
			int qloc, idx, len;
			bool isLess;
			CandidateLoc(int qloc, int idx, int len, bool isLess) : qloc(qloc), idx(idx), len(len), isLess(isLess) {}

			bool operator<(const CandidateLoc& a) const {
				return len < a.len;
			}
		};


		int64_t get_memory_usage()
		{
			int64_t ret = sizeof(*this);
        	ret += sizeof(uint32_t) * int64_t(nPnts);        //checked
			ret += sizeof(std::vector<DataIdx>) * sorted_idx.capacity() + sizeof(DataIdx) * int64_t(nPnts) * dim;  //sorted_idx
			ret += sizeof(std::vector<int32_t>) * next_link.capacity() + sizeof(int32_t) * int64_t(nPnts) * dim;  //next
			return ret;
		}

		template<typename F>
		void for_candidates(int nCandidates, const std::vector<int32_t>& query, const F& f) 
		{
			const int32_t* queryp = (const int32_t*)&query[0];

			std::vector<CandidateLoc> pool;
			pool.reserve(nSearchLoc * 2 + nCandidates*nSearchLoc);
			std::priority_queue<CandidateLoc, std::vector<CandidateLoc>> candidates(std::less<CandidateLoc>() , std::move(pool));
			
            // inque.clear();
            checked.clear();
            int checkCounter = 0;
            //call the callback function if not checked and increase the counter
            const auto& tryCheck = [&](int dataidx){
				if (!checked.isMarked(dataidx)) {
                    checked.mark(dataidx);
					f(dataidx);
                    checkCounter++;
				}
            };

			//it could be possible to maintain only k candidates
			const auto& addCandidates = [&](int qloc, int idx, int plen, bool isLess) {
				// printf("    %d, %d, %d, %d\n", qloc, idx, plen, isLess);
				// int qlocidx = idx*dim + qloc;
				// if (!inque.isMarked(qlocidx)) {
					// printf("%d, %d, %d, %d\n", qloc, idx, plen, isLess);
				candidates.emplace(qloc, idx, plen, isLess);
					// inque.mark(qlocidx);
				// }
			};
			const auto& addCandidates2 = [&](int qloc, int idx) {
				// int qlocidx = idx*dim + qloc;
				// if (!inque.isMarked(qlocidx)) {
				const int32_t* datap = get_datap(qloc, idx);
				auto [plen, isLess] = match_util(queryp, datap, qloc, 0);
				// printf("    addCandidates2   %d, %d, plen=%d, isless=%d\n", qloc, idx, plen, isLess);
				candidates.emplace(qloc, idx, plen, isLess);
					// inque.mark(qlocidx);
				// }
			};

			// printf("------------\n");

			//get the first location of query
			// MyTimer::pusht();
			auto [curidx, lowlen, highlen] = get_loc(queryp, 0);
			addCandidates(0, curidx, lowlen, false);
			addCandidates(0, curidx+1, highlen, true);
			// double t0 =  MyTimer::popt();

			MyTimer::pusht();
			// int potential = 1;
			for(int d_=1;d_<nSearchLoc;d_++){
				int d = d_*step;
				int lowidx = next_link[d-step][curidx];
				int highidx = next_link[d-step][curidx+1];

				// const auto& sorted_idx_d = sorted_idx[d];
				// printf("  low, lowlen, high, highlen = %d, %d, %d, %d\n", lowidx, lowlen, highidx, highlen);

				if(lowlen < step || highlen < step){
					// potential = nCandidates;
					std::tie(curidx, lowlen, highlen) = get_loc(queryp, d);
				} else{
					if(lowlen!=dim){
						lowlen -= step;
					}
					if(highlen!=dim){
						highlen -= step;
					}
					// potential = potential*(highidx-lowidx)+1;
					// std::tie(curidx, lowlen, highlen) = get_loc_scan(queryp, d, lowidx, lowlen, highidx, highlen);
					std::tie(curidx, lowlen, highlen) = get_loc_mixed(queryp, d, lowidx, lowlen, highidx, highlen);
				}
				// printf("  curidx, lowlen, highlen=%d, %d, %d\n", curidx, lowlen, highlen);
				// printf("  push %d, %d, %d, %d\n", d, curidx, lowlen, false);
				// printf("  push %d, %d, %d, %d\n", d, curidx+1, highlen, true);
				// addCandidates(d, curidx, lowlen, false);
				// addCandidates(d, curidx+1, highlen, true);
				// if(potential>=nCandidates){
				// 	potential=1;
				for(int i=curidx;i>=0 && curidx-i<nCandidates;--i){
					int matchingIdx = getidx(d, i);
					tryCheck(matchingIdx);
				}
				for(int i=curidx+1;i<nPnts && i-curidx-1<nCandidates;i++){
					int matchingIdx = getidx(d, i);
					tryCheck(matchingIdx);
				}
				// }
			}
		}
	};


	class LCCS_LSH
	{
	public:
		LCCS_LSH(int M, int step)
			:M(M), step(step), lccs(M*step, step)
		{
			
		}

		int M;
		int step;
		int nPnts;

		// NDArray<2, int32_t> combined_codes;
		// HashCombinator<int32_t> hc;

        void build(NDArray<2, int32_t> &codes)
        {
			nPnts = codes.lens[0];

			lccs.build(codes);
        }
		LCS_SORT_INT lccs;


		int64_t get_memory_usage()
		{
			return lccs.get_memory_usage();
		}

		template<typename F>
		void for_candidates(int nCandidates, const std::vector<int32_t>& query, const F& f) 
		{
			lccs.for_candidates(nCandidates, query, f);
		}
	};

}
