#include "lcsb.h"

using namespace std;

// int MATCH_CNT = 0;

uint8_t clz16(						// count leading zeros (clz)	
	uint16_t x)							// input uint16_t value
{
	uint8_t i = 0;
	if (x == 0) {
		i = 16;
	}
	else {
		for (i = 0; ((x & 0x8000) == 0); i++, x <<= 1);
	}
	return i;
}

float uniform(						// r.v. from Uniform(min, max)
	float min,							// min value
	float max)							// max value
{
	int   num  = rand();
	float base = (float) RAND_MAX - 1.0F;
	float frac = ((float) num) / base;

	return (max - min) * frac + min;
}

float gaussian(						// r.v. from Gaussian(mean, sigma)
	float mean,							// mean value
	float sigma)						// std value
{
	float v1 = -1.0f;
    float v2 = -1.0f;
	float s  = -1.0f;
	float x  = -1.0f;

	do {
		v1 = 2.0F * uniform(0.0F, 1.0F) - 1.0F;
		v2 = 2.0F * uniform(0.0F, 1.0F) - 1.0F;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.0F);
	x = v1 * sqrt (-2.0F * log (s) / s);

	x = x * sigma + mean; 			// x is distributed from N(0, 1)
	return x;
}

// -----------------------------------------------------------------------------
int cmp_item(						// item compare func for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2) 					// 2nd element
{
	Item *x = (Item*) e1;
	Item *y = (Item*) e2;

	int n     = (int) ceil(x->n_ / 64.0f);
	int shift = (int) floor(x->pos_ / 64.0f);
	int start = x->pos_ % 64;

	uint64_t x_i = x->key_[shift] << start;
	uint64_t y_i = y->key_[shift] << start;	

	if (x_i < y_i) return -1;
	else if (x_i > y_i) return 1;
	else {
		for (int i = 1; i < n; ++i) {
			int idx = (shift + i) % n;
			if (x->key_[idx] < y->key_[idx]) return -1;
			else if (x->key_[idx] > y->key_[idx]) return 1;
		}

		if (start > 0) {
			// x_i = (x->key_[shift] >> (64 - start)) << (64 - start);
			// y_i = (y->key_[shift] >> (64 - start)) << (64 - start);
			x_i = x->key_[shift] >> (64 - start);
			y_i = y->key_[shift] >> (64 - start);

			if (x_i < y_i) return -1;
			else if (x_i > y_i) return 1;
		}

		if (x->id_ < y->id_) return -1;
		else if (x->id_ > y->id_) return 1;
		else return 0;
	}
}

// -----------------------------------------------------------------------------
//  Longest Common Substring Bucketing (or simply LCSB) is used to solve the 
//  problem of Maximum Cosine Similarity Search (MCSS)
// -----------------------------------------------------------------------------
LCSB::LCSB(							// constructor
	int n,								// cardinality of dataset
	int d,								// dimensionality of dataset
	int m,								// number of hash tables
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	m_     = m;
	data_  = data;

	// -------------------------------------------------------------------------
	//  generate random projection vectors
	// -------------------------------------------------------------------------
	gen_random_vectors();

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	bulkload();
}

// -----------------------------------------------------------------------------
LCSB::~LCSB()						// destructor
{
	if (proj_ != NULL) {
		for (int i = 0; i < m_; ++i) {
			delete[] proj_[i];	proj_[i] = NULL;
		}
		delete[] proj_;	proj_ = NULL;
	}

	if (list_ != NULL) {
		for (int i = 0; i < m_; ++i) {
			delete[] list_[i]; list_[i] = NULL;
		}
		delete[] list_; list_ = NULL;
	}

	if (hash_key_ != NULL) {
		for (int i = 0; i < n_pts_; ++i) {
			delete[] hash_key_[i]; hash_key_[i] = NULL;
		}
		delete[] hash_key_; hash_key_ = NULL;
	}

	if (table16_ != NULL) {
		delete[] table16_; table16_ = NULL;
	}

	if (checked_ != NULL) {
		delete[] checked_; checked_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void LCSB::gen_random_vectors()		// generate random projection vectors
{
	proj_ = new float*[m_];
	for (int i = 0; i < m_; ++i) {
		proj_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			proj_[i][j] = gaussian(0.0f, 1.0f);
		}
	}
}

// -----------------------------------------------------------------------------
void LCSB::bulkload()				// bulkloading
{
	// -------------------------------------------------------------------------
	//  initialize lookup table for all uint16_t values and checked array
	// -------------------------------------------------------------------------
	int size = 1 << 16;
	table16_ = new uint8_t[size];
	for (int i = 0; i < size; ++i) {
		table16_[i] = clz16(i);
		// printf("i = %5d, val = %d\n", i+1, table16_[i]);
	}

	cnt_     = 1;
	checked_ = new uint32_t[n_pts_];
	memset(checked_, 0, sizeof(uint32_t) * n_pts_);
	printf("Finish lookup table and checked array\n");

	// -------------------------------------------------------------------------
	//  calculate hash code after random projection
	// -------------------------------------------------------------------------
	bool *hash_code = new bool[m_];
	int n = (int) ceil(m_ / 64.0f);

	hash_key_ = new uint64_t*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < m_; ++j) {
			hash_code[j] = calc_hash_code(j, data_[i]);
		}
		hash_key_[i] = compress_hash_code(n, hash_code);
	}
	delete[] hash_code; hash_code = NULL;
	printf("Finish random projection\n");


	//init extents
	init_extent();

	// -------------------------------------------------------------------------
	//  build and sort the bucket list
	// -------------------------------------------------------------------------
	list_ = new Item*[m_];
	for (int i = 0; i < m_; ++i) {
		list_[i] = new Item[n_pts_];
		for (int j = 0; j < n_pts_; ++j) {
			list_[i][j].key_ = hash_key_[j];
			list_[i][j].n_   = m_;
			list_[i][j].id_  = j;
			list_[i][j].pos_ = i;		
		}
		qsort(list_[i], n_pts_, sizeof(Item), cmp_item);
	}
	printf("Finish bucket list\n\n");

	// -------------------------------------------------------------------------
	//  display the bucket list after sorting
	// -------------------------------------------------------------------------
	// for (int i = 0; i < m_; ++i) {
	// 	for (int j = 0; j < n_pts_; ++j) {
	// 		printf("[ ");
	// 		for (int u = 0; u < list_[i][j].n_; ++u) {
	// 			printf("%llu ", list_[i][j].key_[u]);
	// 		}
	// 		printf("], id = %d\n", list_[i][j].id_);
	// 	}
	// 	printf("\n");
	// 	break;
	// }
}

// -----------------------------------------------------------------------------
bool LCSB::calc_hash_code(			// calc hash code after random projection
	int   id,							// projection vector id
	const float *data)					// input data
{
	return calc_inner_product(dim_, proj_[id], data) >= 0 ? true : false;
}

// -----------------------------------------------------------------------------
uint64_t* LCSB::compress_hash_code(	// compress hash code with 64 bits
	int  n,								// compress length
	bool *hash_code)					// input hash code
{
	uint64_t *hash_key = new uint64_t[n];
	int shift = 0;
	for (int i = 0; i < n; ++i) {
		uint64_t val = 0;
		int size = (i == n-1 && m_%64 != 0) ? (m_ % 64) : 64;

		for (int j = 0; j < size; ++j) {
			int idx = (j + shift) % m_;
			if (hash_code[idx]) val = val | ((uint64_t) 1 << (63-j));
		}
		hash_key[i] = val;
		shift = (shift + size) % m_;
	}
	return hash_key;
}

// -----------------------------------------------------------------------------
int LCSB::mcss(						// mcss
	int   top_k,						// top-k number of data objects
	int   check_k,						// check-k number of data objects
	const float *query,					// input query
	MaxK_List *list)					// top-k MC results (return)
{
	// -------------------------------------------------------------------------
	//  calculate the hash code and get the compressed hash key of query
	// -------------------------------------------------------------------------
	bool *hash_code_q = new bool[m_];
	for (int i = 0; i < m_; ++i) {
		hash_code_q[i] = calc_hash_code(i, query);
	}

	int n = (int) ceil(m_ / 64.0f);
	uint64_t *hash_key_q = compress_hash_code(n, hash_code_q);

	// -------------------------------------------------------------------------
	//  initialize the priority queue for candidate generation
	// -------------------------------------------------------------------------
	priority_queue<Elem> cand;
	int id    = -1;
	int loc   = -1;
	int left  = -1;
	int right = -1;
	int match = -1;

	double t1 = MyTimer::measure([&](){
	int step = (int) ceil(m_ * 0.05f);
	// int step = 1;
	for (loc = 0; loc < m_; loc = loc + step) {
		int pos = binary_search(loc, hash_key_q, match);
		// printf("i = %3d, pos = %5d, match = %d\n", loc, pos, match);

		if (match > 0) {
			Elem obj(match, loc, pos-1, pos+1, list_[loc][pos].id_);
			cand.push(obj);
		}
	}
	});

	// -------------------------------------------------------------------------
	//  m-way merge sort and verification
	// -------------------------------------------------------------------------
	int count       = 0;
	int left_match  = -1;
	int right_match = -1;
	int cand_size   = check_k + top_k;
	cand_size = cand_size > n_pts_ ? n_pts_ : cand_size;

	printf("----------\n");

	// printf("before:que.size()=%d\n", cand.size());
	double t2 = MyTimer::measure([&](){
	while (!cand.empty()) {
		Elem top = cand.top();
		id    = top.id_;
		loc   = top.loc_;
		left  = top.left_;
		right = top.right_;
		match = top.key_;

		cand.pop();

		printf("    %d, %d, %d, %d, %d\n", id, loc, left, right, match);

		// ---------------------------------------------------------------------
		//  verification
		// ---------------------------------------------------------------------
		if (checked_[id] < cnt_) {
			list->insert(calc_cosangle(dim_, data_[id], query), id + 1);
			checked_[id] = cnt_;
			if (++count >= cand_size) break;
		}

		// ---------------------------------------------------------------------
		//  get the next 2 possible candidates
		// ---------------------------------------------------------------------
		left_match = -1;			// case 1: left position
		if (left >= 0) cmp_hash_key(loc, left, hash_key_q, left_match);
		if (left_match == match) {
			Elem obj(left_match, loc, left-1, right, list_[loc][left].id_);
			cand.push(obj);
			continue;
		}

		right_match = -1;			// case 2: right position
		if (right < n_pts_) cmp_hash_key(loc, right, hash_key_q, right_match);
		if (right_match > left_match) {
			Elem obj(right_match, loc, left, right+1, list_[loc][right].id_);
			cand.push(obj);
		}
		else {
			Elem obj(left_match, loc, left-1, right, list_[loc][left].id_);
			cand.push(obj);
		}
	}
	});
	// printf("    %f, %f\n", t1, t2);
	// printf("match_cnt=%d\n", MATCH_CNT);
	// printf("\n");

	// -------------------------------------------------------------------------
	//  update <cnt_> and <checked_>
	// -------------------------------------------------------------------------
	cnt_++;
	if (cnt_ == 4294967295U) {
		cnt_ = 1;
		memset(checked_, 0, sizeof(uint32_t) * n_pts_);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	while (!cand.empty()) cand.pop();
	delete[] hash_key_q;
	delete[] hash_code_q; 
	
	return 0;
}



//O(dim * 2^K + dim * n), which is accectable if K is chosen carefully 
void LCSB::init_extent()
{
	extents.resize(m_);
	for (int i = 0; i < extents.size(); i++) {
		std::fill(extents[i].begin(), extents[i].end(), 0);
	}

	for (int i = 0; i < n_pts_; i++) {
		const uint64_t* key_i = hash_key_[i];
		for (int d = 0; d < m_; d++) {
			//K bits prefix at location-d
			uint64_t prefix_datai = prefix(key_i, d, K);
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
	uint64_t prefixQ = prefix(hash_key, loc);
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

	Item &x = list_[loc][pos];
	match = 0;

	int n     = (int) ceil(x.n_ / 64.0f);
	int shift = (int) floor(x.pos_ / 64.0f);
	int last  = (64 - x.n_ % 64) % 64;
	int start = x.pos_ % 64;
	
	uint64_t x_i = x.key_[shift] << start;
	uint64_t y_i = hash_key[shift] << start;
	match += table_lookup(x_i ^ y_i);

	if (x_i < y_i) return -1;
	else if (x_i > y_i) return 1;
	else {
		if (shift == n-1) match -= (last + start);
		else match -= start;

		for (int i = 1; i < n; ++i) {
			int idx = (shift + i) % n;
			match += table_lookup(x.key_[idx] ^ hash_key[idx]);

			if (x.key_[idx] < hash_key[idx]) return -1;
			else if (x.key_[idx] > hash_key[idx]) return 1;
			else if (idx == n-1) match -= last;
		}

		if (start > 0) {
			x_i = (x.key_[shift]   >> (64 - start)) << (64 - start);
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
	if (x == 0) return 64;

	int match = 0;
	uint16_t v = 0;

	v = (x >> 48) & 0xffff;
	match += table16_[v];
	if (v > 0) return match;

	v = (x >> 32) & 0xffff;
	match += table16_[v];
	if (v > 0) return match;

	v = (x >> 16) & 0xffff;
	match += table16_[v];
	if (v > 0) return match;

	v = x & 0xffff;
	match += table16_[v];
	if (v > 0) return match;
}