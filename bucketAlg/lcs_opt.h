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


namespace mylcs 
{//LCS by sort
	class LCS_SORT
	{
	public:
		LCS_SORT(int dim, int step)
            :dim(dim), step(step)
		{
			assert(dim%64==0);
			nSearchLoc = ((dim-1)/step)+1;
		}

        void build(NDArray<2, uint64_t> &data)
        {
            dim = dim;
            // datap = &data;
            datap = data.to_ptr();
			// nPnts = datap->size();
			nPnts = data.lens[0];

            checked.resize(nPnts);

			initExtent();
			initSortedIdx();
			initNextLink();
        }

		struct DataIdx
		{
			static const int StartLoc = 32;
			static const int EndLoc = StartLoc + 32;

			uint32_t prefix_cache;
			uint32_t idx;

			DataIdx(): prefix_cache(0), idx(0){}
			DataIdx(uint32_t p, uint32_t i): prefix_cache(p), idx(i){}
		};

		int dim;
		int nPnts;
		// double alpha;
		int step;
		int nSearchLoc;
		// std::vector<std::vector<uint64_t> > *datap;
		uint64_t** datap;
		// std::vector<std::vector<int32_t> > sorted_idx;
		std::vector<std::vector<DataIdx> > sorted_idx;
		std::vector<std::vector<int32_t> > next_link;

		// static const int TableSize = (1 << K);
		// std::vector<std::array<int32_t, TableSize + 1> > extents;
		std::vector<std::vector<int32_t> > extents;
		int logn;

        CountMarkerU<uint16_t> checked;
        // CountMarkerU<uint16_t> inque;

		inline const uint64_t* getDatap(int qloc, int idx)
		{
			// return (const uint64_t*)&datap->at(sorted_idx[qloc][idx])[0];
			// return (const uint64_t*)&datap->at(sorted_idx[qloc][idx].idx)[0];
			return datap[sorted_idx[qloc][idx].idx];
		}
		inline const uint64_t* getDatap(const DataIdx& idx)
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

		inline uint64_t getSigL(const uint64_t* sigs, int loc)
		{
			int loc64 = loc / 64;
			int locmod64 = loc % 64;
			uint64_t r0 = sigs[loc64] << locmod64;
			return r0;
		}
		inline uint64_t getSigR(const uint64_t* sigs, int loc)
		{
			if(loc%64==0){
				return 0;
			}
			int loc64 = loc / 64;
			int locmod64 = loc % 64;
			int loc64p = (loc64 + 1) % dim64();
			uint64_t r1 = sigs[loc64p] >> (64 - locmod64);
			return r1;
		}
		//return 64 bits integer sigs[loc:loc+64]
		inline uint64_t getSig(const uint64_t* sigs, int loc)
		{
			uint64_t r0 = getSigL(sigs, loc);
            uint64_t r1 = getSigR(sigs, loc);
			return r0 ^ r1;
		}

		//return 64 bits integer reverse(sigs[loc-64:loc])
		//inline uint64_t getBSig()
		inline uint64_t prefix(uint64_t u, int n)
		{
			return u >> (64 - n);
		}

		inline int prefix(const uint64_t* sigs, int loc, int n)
		{
			return prefix(getSig(sigs, loc), n);
		}

		//O(dim * 2^K + dim * n), which is accectable if K is chosen carefully 
		void initExtent()
		{
			logn = ceil(log(nPnts)/log(2.))-2;
			int TableSize = (1<<logn)+1;
			extents.resize(dim);
			for (int i = 0; i < extents.size(); i++) {
				extents[i].resize(TableSize);
				std::fill(extents[i].begin(), extents[i].end(), 0);
			}

			for (int i = 0; i < nPnts; i++) {
				// auto dataip = (const uint64_t*)&datap->at(i)[0];
				auto dataip = datap[i];
				for (int d = 0; d < dim; d++) {
					//K bits prefix at location-d
					uint64_t prefix_datai = prefix(dataip, d, logn);
					extents[d][prefix_datai + 1]++;
				}
			}

			//sum the extent
			for (int d = 0; d < dim; d++) {
				for (int i = 0; i < TableSize; i++) {
					extents[d][i + 1] += extents[d][i];
				}
			}
		}

		//bucket-sort leveraging extents
		void initSortedIdx()
		{
			int TableSize = (1<<logn)+1;
			sorted_idx.resize(dim);
			for (int d = 0; d < dim; d++) {
				sorted_idx[d].resize(nPnts);

				std::vector<int32_t> extent_copy = extents[d];
				for (int i = 0; i < nPnts; i++) {
					// auto dataip = (const uint64_t*)&datap->at(i)[0];
					auto dataip = datap[i];
					uint64_t sig = getSig(dataip, d);
					uint64_t p = prefix(sig, logn);
					sorted_idx[d][extent_copy[p]++] = {uint32_t(sig>>DataIdx::StartLoc), i};
				}

				// const auto& cmpf = [d, this](int a, int b) -> bool {
				// 	const auto &da = &(datap->at(a)[0]);
				// 	const auto &db = &(datap->at(b)[0]);
				// 	// return cmpPrefix(da, db, d);
				// 	auto [len, isLess] = matchUtil(da, db, d, 0); 
				// 	return isLess;
				// };
				const auto& cmpf = [d, this](const DataIdx& a, const DataIdx& b) -> bool {
					auto [len, isLess] = matchUtil(a, b, d, 0); 
					return isLess;
				};

				for (int i = 0; i < TableSize; i++) {
					int l = extents[d][i];
					int h = extents[d][i + 1];
					sort(sorted_idx[d].begin() + l, sorted_idx[d].begin() + h, cmpf);
				}
			}
		}

		void initNextLink()
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
		std::tuple<int, int, int> getLocScan(const uint64_t* q, int loc, int low, int lowlen, int high, int highlen)
		{
			int lastlen = lowlen;
			int minlen = std::min(lowlen, highlen);
			for(int i=low+1;i<high;i++){
				// const uint64_t* dpi = (const uint64_t*)&(datap->at(sorted_idx_loc[i])[0]);
				const uint64_t* dpi = getDatap(loc, i);
				auto [ilen, iless] = matchUtil(q, dpi, loc, minlen);

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
		std::tuple<int, int, int> getLocMixed(const uint64_t* query, int qloc, int low, int lowlen, int high, int highlen)
		{
			static const int SCAN_SIZE = 8;
			// const auto& sorted_idx_qloc = sorted_idx[qloc];

			//return datap->at(idx[qloc][a]) < query
			const auto& cmp = [&](int a, int match=0) -> std::tuple<int, bool> {
				// const uint64_t* dpa = (const uint64_t*)&(datap->at(sorted_idx_qloc[a])[0]);
				const uint64_t* dpa = getDatap(qloc, a);
				return matchUtil(query, dpa, qloc, match);
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

			return getLocScan(query, qloc, low, lowlen, high, highlen);
		}

		//return  (lowidx, lowlen, highlen) that q \in [lowidx, highidx) or lowidx==0 or highidx==N-1
		std::tuple<int, int, int> getLoc(const uint64_t* query, int qloc)
		{
			uint64_t prefixq = prefix(query, qloc, logn);
			int low = extents[qloc][prefixq];
			int high = extents[qloc][prefixq + 1]-1;
			if(low ==high){
				if(low!=nPnts-1){
					high = low+1;
				} else{
					low = high-1;
				}
			}
			// const auto& sorted_idx_qloc = sorted_idx[qloc];
			// const uint64_t* dplow = (const uint64_t*)&(datap->at(sorted_idx_qloc[low])[0]);
			// const uint64_t* dphigh = (const uint64_t*)&(datap->at(sorted_idx_qloc[high])[0]);
			const uint64_t* dplow = getDatap(qloc, low);
			const uint64_t* dphigh = getDatap(qloc, high);
			auto [lowlen, lowisless] = matchUtil(query, dplow, qloc, 0);
			auto [highlen, highisless] = matchUtil(query, dphigh, qloc, 0);

			return getLocMixed(query, qloc, low, lowlen, high, highlen);

		}


		//matchUtil :: [sigs] a-> [sigs] b->int loc-> (int, bool)
		//return matchlen and bool(x<y)
		std::tuple<int, bool> matchUtil(const uint64_t* x, const uint64_t* y, int loc, int match)
		{
			//to make sure the order is correct
			if(match>0){
				--match;
			}
			int matchloc = (loc+match)%dim;

			uint64_t x_i = getSigL(x, matchloc);
			uint64_t y_i = getSigL(y, matchloc);
			match += get_num_prefix(~x_i^y_i);

			if (x_i != y_i) return std::make_tuple(match, x_i < y_i);
			match -= matchloc%64;

			int restlen = dim - match;
			for (int i = 0; i*64 < restlen; i++) {
				int idx = (matchloc/64 + 1 + i) % dim64();
				x_i = x[idx];
				y_i = y[idx];
				match += get_num_prefix(~x_i^y_i);

				if (x_i!=y_i) return std::make_tuple(match, x_i < y_i);
			}
			//equal
			match = dim;
			return std::make_tuple(match, false);
		} 

		std::tuple<int, bool> matchUtil(const DataIdx& x, const DataIdx& y, int loc, int match)
		{
			//to make sure the order is correct
			if(match>=DataIdx::StartLoc && x.prefix_cache != y.prefix_cache){
				int retlen = DataIdx::StartLoc + get_num_prefix(uint32_t(~x.prefix_cache^y.prefix_cache));
				return std::make_tuple(retlen, x.prefix_cache < y.prefix_cache);
			}
			// const uint64_t* xp = &datap->at(x.idx)[0];
			// const uint64_t* yp = &datap->at(y.idx)[0];
			const uint64_t* xp = datap[x.idx];
			const uint64_t* yp = datap[y.idx];
			return matchUtil(xp, yp, loc, match);
		} 

		struct CandidateLoc {
			int qloc, idx, len;
			bool isLess;
			CandidateLoc(int qloc, int idx, int len, bool isLess) : qloc(qloc), idx(idx), len(len), isLess(isLess) {}

			bool operator<(const CandidateLoc& a) const {
				return len < a.len;
			}
		};



		template<typename F>
		void for_candidates(int nCandidates, const std::vector<uint64_t>& query, const F& f) 
		{
			const uint64_t* queryp = (const uint64_t*)&query[0];

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
				const uint64_t* datap = getDatap(qloc, idx);
				auto [plen, isLess] = matchUtil(queryp, datap, qloc, 0);
				// printf("    addCandidates2   %d, %d, plen=%d, isless=%d\n", qloc, idx, plen, isLess);
				candidates.emplace(qloc, idx, plen, isLess);
					// inque.mark(qlocidx);
				// }
			};

			// printf("------------\n");

			//get the first location of query
			// MyTimer::pusht();
			auto [curidx, lowlen, highlen] = getLoc(queryp, 0);
			addCandidates(0, curidx, lowlen, false);
			addCandidates(0, curidx+1, highlen, true);
			// double t0 =  MyTimer::popt();

			MyTimer::pusht();
			for(int d_=1;d_<nSearchLoc;d_++){
				int d = d_*step;
				int lowidx = next_link[d-step][curidx];
				int highidx = next_link[d-step][curidx+1];

				// const auto& sorted_idx_d = sorted_idx[d];
				// printf("  low, lowlen, high, highlen = %d, %d, %d, %d\n", lowidx, lowlen, highidx, highlen);

				if(lowlen < step || highlen < step){
					std::tie(curidx, lowlen, highlen) = getLoc(queryp, d);
				} else{
					if(lowlen!=dim){
						lowlen -= step;
					}
					if(highlen!=dim){
						highlen -= step;
					}
					// std::tie(curidx, lowlen, highlen) = getLocScan(queryp, d, lowidx, lowlen, highidx, highlen);
					std::tie(curidx, lowlen, highlen) = getLocMixed(queryp, d, lowidx, lowlen, highidx, highlen);
				}
				// printf("  push %d, %d, %d, %d\n", d, curidx, lowlen, false);
				// printf("  push %d, %d, %d, %d\n", d, curidx+1, highlen, true);
				addCandidates(d, curidx, lowlen, false);
				addCandidates(d, curidx+1, highlen, true);
			}
			double t1 = MyTimer::popt();

			// printf("~~before  que.size()=%d\n", candidates.size());
			//merge
			// t2= MyTimer::measure([&](){
			MyTimer::pusht();
			while (checkCounter<nCandidates && candidates.size()>0) {
				auto cl = candidates.top();
				candidates.pop();				
				// CandidateLoc* clp = candidates.top();
				// candidates.pop();
				// auto& cl = *clp;

				// int matchingIdx = sorted_idx[cl.qloc][cl.idx];
				int matchingIdx = getidx(cl.qloc, cl.idx);
				tryCheck(matchingIdx);
				printf("           %d: loc=%d, idx=%d, len=%d, isless=%d\n", matchingIdx, cl.qloc, cl.idx, cl.len, cl.isLess);

				if (cl.idx > 0 && !cl.isLess) {
					addCandidates2(cl.qloc, cl.idx - 1);
				}
				if (cl.idx < nPnts - 1 && cl.isLess) {
					addCandidates2(cl.qloc, cl.idx + 1);
				}
			}
			double t2 = MyTimer::popt();

			printf("    t1, t2:%f, %f\n", t1, t2);
			// printf("    t10, t11, t12:%f, %f, %f\n", t10, t11, t12);
			// printf("~~after  que.size()=%d\n", candidates.size());
		}
	};

}
