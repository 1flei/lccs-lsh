#include "sort_lcs.h"

using namespace std;

SALCS::SALCS(int dim, std::vector<std::vector<SigType> > &data)
	:dim(dim), datap(&data)
{
	idx.resize(dim);
	for (int i = 0; i < dim; i++) {
		idx[i].resize(data.size());
		for (int j = 0; j < idx[i].size(); j++) {
			idx[i][j] = j;
		}
	}

	for (int i = 0; i < dim; i++) {
		const auto& cmp = [&](int a, int b) -> bool {
			const vector<SigType> &da = data[a];
			const vector<SigType> &db = data[b];
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

		sort(idx[i].begin(), idx[i].end(), cmp);
	}
}


int SALCS::getLoc(std::vector<SigType>& query, int qloc)
{
	//return datap->at(idx[qloc][a]) <= query
	const auto& cmp = [&](int a) -> bool {
		const vector<SigType> &da = datap->at(idx[qloc][a]);
		for (int j = 0; j < da.size(); j++) {
			int loc = (j + qloc) % dim;
			if (da[loc] < query[loc]) {
				return true;
			}
			else if (da[loc] > query[loc]) {
				return false;
			}
		}
		return true;
	};

	//make sure that low <=query < high
	int low = 0;
	int high = idx[qloc].size();

	while (low < high -1) {
		//binary search
		int mid = (low + high) / 2;

		//mid <= query, so query in [mid, high)
		if (cmp(mid)) {
			low = mid;
		}
		else {		//otherwise query in [low, mid)
			high = mid;
		}
	}

	//now low should be the location
	return low;
}



int SALCS::matchLen(std::vector<SigType>& a, std::vector<SigType>& b, int qloc)
{
	for (int i = 0; i < dim; i++){
		int loc = (i + qloc) % dim;
		if (a[loc] != b[loc]) {
			return i;
		}
	}
	return dim;
}