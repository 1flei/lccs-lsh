#include "headers.h"

// -----------------------------------------------------------------------------
KL_SRP::KL_SRP(						// constructor
	int n,								// cardinality of dataset
	int d,								// dimensionality of dataset
	int k,								// number of projections per table
	int l,                              // number of hash tables
    const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	k_     = k;
	l_     = l;	
	m_     = (int) ceil(k_ / 31.0f);
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
KL_SRP::~KL_SRP()					// destructor
{
	if (proj_ != NULL) {
		for (int i = 0; i < l_; ++i) {
            for (int j = 0; j < k_; ++j) {
                delete[] proj_[i][j]; proj_[i][j] = NULL;
            }
			delete[] proj_[i];	proj_[i] = NULL;
		}
		delete[] proj_;	proj_ = NULL;
	}

	if (standard_hash_ != NULL) {
		delete[] standard_hash_; standard_hash_ = NULL;
	}

	if (!tables_.empty()) {
        tables_.clear(); tables_.shrink_to_fit();
	}

	if (checked_ != NULL) {
		delete[] checked_; checked_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void KL_SRP::gen_random_vectors()	// generate random projection vectors
{
	// -------------------------------------------------------------------------
	//  generate random projection vectors
	// -------------------------------------------------------------------------
	proj_ = new float**[l_];
	for (int i = 0; i < l_; ++i) {
		proj_[i] = new float*[k_];
		for (int j = 0; j < k_; ++j) {
			proj_[i][j] = new float[dim_];
            for (int u = 0; u < dim_; ++u) {
                proj_[i][j][u] = gaussian(0.0f, 1.0f);
            }
		}
	}

	// -------------------------------------------------------------------------
	//  generate standard hash function
	// -------------------------------------------------------------------------
	uint32_t max_hash_rnd = 536870912U; // 2^29

	standard_hash_ = new uint32_t[m_];
	for (int i = 0; i < m_; ++i) {
		standard_hash_[i] = uniform32(1, max_hash_rnd);
	}
}

// -----------------------------------------------------------------------------
void KL_SRP::bulkload()				// bulkloading
{
	// -------------------------------------------------------------------------
	//  initialize counter and checked array
	// -------------------------------------------------------------------------
	cnt_     = 1;
	checked_ = new uint32_t[n_pts_];
	memset(checked_, 0, sizeof(uint32_t) * n_pts_);

	// -------------------------------------------------------------------------
	//  allocate the space for hash tables
	// -------------------------------------------------------------------------
	for (int i = 0; i < l_; ++i) {
		unordered_map<int, vector<int> > table;
		tables_.push_back(table);
	}

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < l_; ++j) {
			int val = calc_hash_value(j, data_[i]);
			tables_[j][val].push_back(i);
		}
	}
}

// -----------------------------------------------------------------------------
int KL_SRP::calc_hash_value(		// calculate hash value of input data
	int   hash_id,                  	// hash table id
	const float *data)					// input data
{
	uint32_t largest32 = 4294967295U; // 2^32 - 1
	uint32_t prime = 4294967291U;	// 2^32 - 5

	uint64_t val = 0;
	int shift = 0;
	for (int i = 0; i < m_; ++i) {
		// ---------------------------------------------------------------------
		//  calculate the hash key for every 31 bits
		// ---------------------------------------------------------------------
		int hash_key = 0;
		int size = (i == m_-1 && k_%31 != 0) ? (k_ % 31) : 31;

		for (int j = 0; j < size; ++j) {
			int   idx = (j + shift) % k_;
			float ip  = calc_inner_product(dim_, proj_[hash_id][idx], data);
			if (ip >= 0) hash_key = hash_key | (1 << j);
		}
		shift = (shift + size) % k_;

		// ---------------------------------------------------------------------
		//  calculate the hash value with incremental update
		// ---------------------------------------------------------------------
		val = val + (uint64_t) hash_key * (uint64_t) standard_hash_[i];
		val = (val & largest32) + 5 * (val >> 32);
		if (val >= prime) {			// fast compute "mod" function
			val = val - prime;
		}
	}
	return (int) val;
}

// -----------------------------------------------------------------------------
int KL_SRP::mcss(					// c-k-AMC search
	int   top_k,						// top-k value
	int   check_k,						// check-k number of data objects
	const float *query,					// input query
	MaxK_List *list)					// top-k MC results  (return)
{
	// -------------------------------------------------------------------------
	//  find candidates which collide with the query in the hash buckets
	// -------------------------------------------------------------------------
	int count = 0;
	int cand_size = check_k + top_k;
	cand_size = cand_size > n_pts_ ? n_pts_ : cand_size;

 	for (int i = 0; i < l_; ++i) {
		int q_val = calc_hash_value(i, query);

		vector<int> &v = tables_[i][q_val];
		int size = v.size();
		for (int j = 0; j < size; ++j) {
			// -----------------------------------------------------------------
			//  verification
			// -----------------------------------------------------------------
			int id = v[j];
			if (checked_[id] < cnt_) {
				list->insert(calc_cosine(dim_, data_[id], query), id + 1);
				checked_[id] = cnt_;
				if (++count >= cand_size) break;
			}
		}
		if (count >= cand_size) break;
	}
	
	// -------------------------------------------------------------------------
	//  update <cnt_> and <checked_>
	// -------------------------------------------------------------------------
	cnt_++;
	if (cnt_ == 4294967295U) {
		cnt_ = 1;
		memset(checked_, 0, sizeof(uint32_t) * n_pts_);
	}

	return 0;
}