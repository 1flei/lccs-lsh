#include "lcs.h"

using namespace std;

int LCSSort::getLoc(const std::vector<SigType>& query, int qloc)
{
	//return datap->at(idx[qloc][a]) <= query
	const auto& cmp = [&](int a) -> bool {
		const vector<SigType> &da = codesp->at(idx[qloc][a]);
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



int LCSSort::matchLen(const std::vector<SigType>& a, const std::vector<SigType>& b, int qloc)
{
	for (int i = 0; i < dim; i++){
		int loc = (i + qloc) % dim;
		if (a[loc] != b[loc]) {
			return i;
		}
	}
	return dim;
}