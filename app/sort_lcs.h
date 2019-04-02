#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <queue>
#include <set>
#include <unordered_set>

typedef int SigType;

//LCS by sort
class SALCS
{
public:
	SALCS(int dim, std::vector<std::vector<SigType> > &data);

	int dim;
	std::vector<std::vector<SigType> > *datap;
	std::vector<std::vector<int> > idx;

	//binary search 
	int getLoc(std::vector<SigType>& query, int qloc);
	int matchLen(std::vector<SigType>& a, std::vector<SigType>& b, int qloc);

	struct CandidateLoc {
		int qloc, idx, len;
		CandidateLoc(int qloc, int idx, int len) : qloc(qloc), idx(idx), len(len) {}

		bool operator<(const CandidateLoc& a) const {
			return len < a.len;
		}
	};

	template<typename F>
	void klcs(std::vector<SigType>& query, int k, const F& f) {
		assert(query.size() == dim);

		std::priority_queue<CandidateLoc> candidates;
		std::set<std::pair<int, int> > inque;

		const auto addToCandidates = [&](int qloc, int idx, int plen) {
			auto p = std::make_pair(qloc, idx);
			if (inque.find(p) == inque.end()) {
				candidates.emplace(qloc, idx, plen);
				inque.emplace(qloc, idx);
			}
		};

		int qloc = 0;
		while (qloc < dim) {
			int idxi = getLoc(query, qloc);
			int prefixLen = matchLen(query, datap->at(idx[qloc][idxi]), qloc);
			if (prefixLen > 0) {
				addToCandidates(qloc, idxi, prefixLen);
				qloc += prefixLen;
			}
			else {
				qloc++;
			}
		}
		//k-way merge for candidates vector

		std::unordered_set<int> checked;

		while (checked.size()<k) {
			auto cl = candidates.top();
			candidates.pop();

			int matchingIdx = idx[cl.qloc][cl.idx];

			if (checked.find(matchingIdx) == checked.end()) {
				checked.insert(matchingIdx);
				f(matchingIdx);
			}
			//push new candidates from cl
			if (cl.len > 1) {
				int idxp1 = getLoc(query, (cl.qloc + 1)%dim);
				addToCandidates((cl.qloc + 1)%dim, idxp1, cl.len-1);
			}
			
			if (cl.idx > 0) {
				int plen = matchLen(query, datap->at(idx[cl.qloc][cl.idx - 1]), cl.qloc);
				addToCandidates(cl.qloc, cl.idx - 1, plen);
			}
			if (cl.idx < datap->size() - 1) {
				int plen = matchLen(query, datap->at(idx[cl.qloc][cl.idx + 1]), cl.qloc);
				addToCandidates(cl.qloc, cl.idx + 1, plen);
			}
		}
	}
};