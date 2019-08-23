#pragma once

#include "../def.h"
#include "../util.h"
#include <vector>
#include <queue>
#include <set>
#include <unordered_set>

class LCSSort
{
public:
    LCSSort(int d) : dim(d) {
        idx.resize(dim);
    };
    int dim;
    const std::vector<std::vector<SigType> > *codesp;
	std::vector<std::vector<int> > idx;
    

	void build(const std::vector<std::vector<SigType> > &codes) {
        codesp = &codes;

        for (int i = 0; i < dim; i++) {
            idx[i].resize(codes.size());
            for (int j = 0; j < idx[i].size(); j++) {
                idx[i][j] = j;
            }
        }

        for (int i = 0; i < dim; i++) {
            const auto& cmp = [&](int a, int b) -> bool {
                const std::vector<SigType> &da = codes[a];
                const std::vector<SigType> &db = codes[b];
                for (int j = 0; j < da.size(); j++) {
                    int loc = (j + i) % dim;
                    if (da[loc] < db[loc]) {
                        return true;
                    }
                    else if (da[loc] > db[loc]) {
                        return false;
                    }
                }
                return false;
            };

            std::sort(idx[i].begin(), idx[i].end(), cmp);
        }
    }

	//binary search 
	int getLoc(const std::vector<SigType>& query, int qloc);
	int matchLen(const std::vector<SigType>& a, const std::vector<SigType>& b, int qloc);

	struct CandidateLoc {
		int qloc, idx, len;
		CandidateLoc(int qloc, int idx, int len) : qloc(qloc), idx(idx), len(len) {}

		bool operator<(const CandidateLoc& a) const {
			return len < a.len;
		}
	};

	template<typename F>
	void for_candidates(int nCandidates, std::vector<SigType>& query, const F& f) {
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
			int prefixLen = matchLen(query, codesp->at(idx[qloc][idxi]), qloc);
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

		while (checked.size()<nCandidates) {
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
				int plen = matchLen(query, codesp->at(idx[cl.qloc][cl.idx - 1]), cl.qloc);
				addToCandidates(cl.qloc, cl.idx - 1, plen);
			}
			if (cl.idx < codesp->size() - 1) {
				int plen = matchLen(query, codesp->at(idx[cl.qloc][cl.idx + 1]), cl.qloc);
				addToCandidates(cl.qloc, cl.idx + 1, plen);
			}
		}
	}
};