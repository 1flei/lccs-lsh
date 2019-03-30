#pragma once
#include "sam.h"
#include <unordered_set>

class LCSIndex
{
public:
	LCSIndex(int dim) : dim(dim), useDelimiter(false) {}
	LCSIndex(int dim, int delimiter) : dim(dim), useDelimiter(true), delimiter(delimiter) {}

	MySAM sam;
	int dim;	//excluding delimiter
	bool useDelimiter;
	int delimiter;

	//assume the x itself can be seperated well 
	void insert(std::vector<int> &x) {
		assert(x.size() == dim);
		sam.extends(x);
		if (useDelimiter) {
			sam.extend(delimiter);
		}
	}
	void insert(const std::initializer_list<int> &cs) {
		std::vector<int> v = cs;
		sam.extends(v);
	}

	template<typename F>
	void klcs(std::vector<int>& query, int k, const F& f) {
		int lenq = query.size();

		std::vector<std::unordered_set<int> > toCheck(lenq+1);

		int curLength = 0;
		const auto fmatch = [&](int idx, int j) {
			curLength++;
			//insert the idx to the corresponding toCheck list
			toCheck[curLength].insert(idx);
		};
		const auto flink = [&](int idx, int j) {
			curLength = sam.len(idx);
			//toCheck[curLength].insert(idx);
		};
		sam.traverse(query, fmatch, flink);

		std::unordered_set<int> checked;

		const auto& fLeaves = [&](int idx)->bool {
			int itemIdx = (sam.len(idx)-1) / dim;
			// printf("idx, len, itemIdx: %d, %d, %d\n", idx, sam.len(idx), itemIdx);
			if (checked.find(itemIdx) == checked.end()) {
				//not checked
				checked.insert(itemIdx);
				f(itemIdx);
			}
			return checked.size() < k;
		};

		//now toCheck has the list to check to get the topk longest common substring
		for (int j = lenq; j >= 0; --j) {
			for (int idxToCheck : toCheck[j]) {
				// printf("match=%d, toCheck=%d\n", j, idxToCheck);
				sam.mapSubtree(fLeaves, idxToCheck);
				if (checked.size() >= k) {
					return;
				}
				int pre = sam.parent(idxToCheck);
				if (pre != -1) {
					toCheck[sam.len(pre)].insert(pre);
				}
			}
		}
	}



	std::vector<int> klcs(std::vector<int>& query, int k) {
		std::vector<int> ret;
		ret.reserve(k);
		const auto& f = [&](int itemIdx) {
			ret.push_back(itemIdx);
		};
		klcs(query, k, f);
		return ret;
	}
};