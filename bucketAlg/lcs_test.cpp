#include "lcs_test.h"

using namespace std;

// int MATCH_CNT = 0;

// -----------------------------------------------------------------------------
int LCSB::cmp_item(						// item compare func for qsort (ascending)
	const Item& e1,						// 1st element
	const Item& e2, 
    int pos) 					// 2nd element
{
	const Item& x = e1;
	const Item& y = e2;

	// int n     = (int) ceil(x.n_ / 64.0f);
	// int shift = (int) floor(x.pos_ / 64.0f);
	// int start = x.pos_ % 64;
	int n     = (int) ceil(m_ / 64.0f);
	int shift = (int) floor(pos / 64.0f);
	int start = pos % 64;

	uint64_t x_i = get_key(x)[shift] << start;
	uint64_t y_i = get_key(y)[shift] << start;	

	if (x_i < y_i) return -1;
	else if (x_i > y_i) return 1;
	else {
		for (int i = 1; i < n; ++i) {
			int idx = (shift + i) % n;
			if (get_key(x)[idx] < get_key(y)[idx]) return -1;
			else if (get_key(x)[idx] > get_key(y)[idx]) return 1;
		}

		if (start > 0) {
			// x_i = (get_key(x)[shift] >> (64 - start)) << (64 - start);
			// y_i = (get_key(y)[shift] >> (64 - start)) << (64 - start);
			x_i = get_key(x)[shift] >> (64 - start);
			y_i = get_key(y)[shift] >> (64 - start);

			if (x_i < y_i) return -1;
			else if (x_i > y_i) return 1;
		}

		if (get_idx(x) < get_idx(y)) return -1;
		else if (get_idx(x) > get_idx(y)) return 1;
		else return 0;
        return 0;
	}
}

// -----------------------------------------------------------------------------
//  Longest Common Substring Bucketing (or simply LCSB) is used to solve the 
//  problem of Maximum Cosine Similarity Search (MCSS)
// -----------------------------------------------------------------------------
LCSB::LCSB(							// constructor
	int n,								// cardinality of dataset
	int m,								// number of hash tables
	int step_)					
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	m_     = m;
    step   = step_;
}

// -----------------------------------------------------------------------------
LCSB::~LCSB()						// destructor
{
	if (list_ != NULL) {
		for (int i = 0; i < m_; ++i) {
			delete[] list_[i]; list_[i] = NULL;
		}
		delete[] list_; list_ = NULL;
	}
	if (checked_ != NULL) {
		delete[] checked_; checked_ = NULL;
	}
}


// void LCSB::build(const std::vector<std::vector<uint64_t> > &data)
void LCSB::build(NDArray<2, uint64_t> &data)
{
    // datap = &data;

    cnt_     = 1;
    checked_ = new uint32_t[n_pts_];
    memset(checked_, 0, sizeof(uint32_t) * n_pts_);
    printf("Finish lookup table and checked array\n");

    // -------------------------------------------------------------------------
    //  calculate hash code after random projection
    // -------------------------------------------------------------------------
    // hash_key_ = new uint64_t*[n_pts_];
    // for (int i = 0; i < n_pts_; ++i) {
    //     hash_key_[i] = (uint64_t*)&(datap->at(i)[0]);
    // }
	hash_key_ = data.to_ptr();


    // -------------------------------------------------------------------------
    //  build and sort the bucket list
    // -------------------------------------------------------------------------
    list_ = new Item*[m_];
    for (int i = 0; i < m_; ++i) {
        list_[i] = new Item[n_pts_];
        for (int j = 0; j < n_pts_; ++j) {
            // list_[i][j].key_ = hash_key_[j];
            // list_[i][j].n_   = m_;
            list_[i][j].id_  = j;
            // list_[i][j].pos_ = i;		
        }
        // qsort(list_[i], n_pts_, sizeof(Item), cmp_item);
        std::sort(list_[i], list_[i]+n_pts_, [&](const Item& e1, const Item& e2){
            return cmp_item(e1, e2, i) < 0;
        });
    }
    printf("Finish bucket list\n\n");

    //init extents
    init_extent();
}

//O(dim * 2^K + dim * n), which is accectable if K is chosen carefully 
void LCSB::init_extent()
{
	logn = ceil(log(n_pts_)/log(2.))-2;
	int TableSize = (1<<logn)+1;
    printf("prefixlen=%d, tablesize=%d\n", logn, TableSize);
	extents.resize(m_);
	for (int i = 0; i < extents.size(); i++) {
		extents[i].resize(TableSize);
		std::fill(extents[i].begin(), extents[i].end(), 0);
	}

	for (int i = 0; i < n_pts_; i++) {
		const uint64_t* key_i = hash_key_[i];
		for (int d = 0; d < m_; d++) {
			//K bits prefix at location-d
			uint64_t prefix_datai = prefix(key_i, d, logn);
			extents[d][prefix_datai + 1]++;
		}
	}

	//sum the extent
	for (int d = 0; d < m_; d++) {
		for (int i = 0; i < TableSize; i++) {
			extents[d][i + 1] += extents[d][i];
		}
	}
}

// -----------------------------------------------------------------------------
int LCSB::binary_search(			// binary search the location of item
	int loc,							// location of the hash code
	uint64_t *hash_key,					// input hash key
	int &match) 						// number of match item (return)
{
	// -------------------------------------------------------------------------
	//  binary search the location
	// -------------------------------------------------------------------------
	uint64_t prefixQ = prefix(hash_key, loc, logn);
	int left   = extents[loc][prefixQ];
	int right  = extents[loc][prefixQ+1];
	// int left   = 0;
	// int right  = n_pts_-1;
	int mid    = 0;
	int result = -1;

	while (left < right) {
		mid = (left + right + 1) / 2;
		result = cmp_hash_key(loc, mid, hash_key, match);

		if (result == 0) return mid;
		else if (result == -1) left = mid;
		else right = mid - 1;
	}
	// assert(left >= 0 && left < n_pts_);

	// -------------------------------------------------------------------------
	//  determine the final position with largest match value (left or right)
	// -------------------------------------------------------------------------
	int left_match = -1;
	cmp_hash_key(loc, left, hash_key, left_match);

	right = left + 1;
	if (right < n_pts_) {
		int right_match = -1;
		cmp_hash_key(loc, right, hash_key, right_match);
		if (right_match > left_match) {
			match = right_match;
			return right;
		}
	}

	match = left_match;
	return left;
}

// -----------------------------------------------------------------------------
int LCSB::cmp_hash_key(				// cmp func for input hash key
	int loc,							// location of hash key
	int pos,							// position of hash key
	uint64_t *hash_key,					// input hash key
	int &match)							// number of match bits (return)
{
	// ++MATCH_CNT;
	// MyTimer tm("cmp");

	Item &x = list_[loc][pos];
	match = 0;

	// int n     = (int) ceil(x.n_ / 64.0f);
	// int shift = (int) floor(x.pos_ / 64.0f);
	// int last  = (64 - x.n_ % 64) % 64;
	// int start = x.pos_ % 64;
	int n     = (int) ceil(m_ / 64.0f);
	int shift = (int) floor(loc / 64.0f);
	int last  = (64 - m_ % 64) % 64;
	int start = loc % 64;
	
	uint64_t x_i = get_key(x)[shift] << start;
	uint64_t y_i = hash_key[shift] << start;
	match += table_lookup(x_i ^ y_i);

	if (x_i < y_i) return -1;
	else if (x_i > y_i) return 1;
	else {
		if (shift == n-1) match -= (last + start);
		else match -= start;

		for (int i = 1; i < n; ++i) {
			int idx = (shift + i) % n;
			match += table_lookup(get_key(x)[idx] ^ hash_key[idx]);

			if (get_key(x)[idx] < hash_key[idx]) return -1;
			else if (get_key(x)[idx] > hash_key[idx]) return 1;
			else if (idx == n-1) match -= last;
		}

		if (start > 0) {
			x_i = (get_key(x)[shift]   >> (64 - start)) << (64 - start);
			y_i = (hash_key[shift] >> (64 - start)) << (64 - start);
			match += table_lookup(x_i ^ y_i);
		
			if (x_i < y_i) return -1;
			else if (x_i > y_i) return 1;
			else match -= (64 - start);
		}
		return 0;
	}
}

// -----------------------------------------------------------------------------
int LCSB::table_lookup(				// table lookup the match value
	uint64_t x)							// input uint64_t value
{
	return get_num_prefix(~x);
}