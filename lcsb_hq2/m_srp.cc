#include "headers.h"

// -----------------------------------------------------------------------------
//  m-Sign-Random Projection LSH (or simply m-srp) is used to solve the 
//  problem of Maximum Cosine Similarity Search (MCSS)
// -----------------------------------------------------------------------------
M_SRP::M_SRP(						// constructor
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
	M_     = (int) ceil(m_ / 64.0f);
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
M_SRP::~M_SRP()						// destructor
{
	if (proj_ != NULL) {
		for (int i = 0; i < m_; ++i) {
			delete[] proj_[i];	proj_[i] = NULL;
		}
		delete[] proj_;	proj_ = NULL;
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
}

// -----------------------------------------------------------------------------
void M_SRP::gen_random_vectors()	// generate random projection vectors
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
void M_SRP::bulkload()				// bulkloading
{
	// -------------------------------------------------------------------------
	//  initialize lookup table for all uint16_t values
	// -------------------------------------------------------------------------
	int size = 1 << 16;
	table16_ = new uint32_t[size];
	for (int i = 0; i < size; ++i) {
		table16_[i] = bit_count(i);
		// printf("%5d: val = %d\n", i, table16_[i]);
	}

	// -------------------------------------------------------------------------
	//  calculate and compress hash code after random projection
	// -------------------------------------------------------------------------
	bool *hash_code = new bool[m_];
	hash_key_ = new uint64_t*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < m_; ++j) {
			hash_code[j] = calc_hash_code(j, data_[i]);
		}
		hash_key_[i] = compress_hash_code(hash_code);
	}
	delete[] hash_code; hash_code = NULL;

	// -------------------------------------------------------------------------
	//  display the bucket list after sorting
	// -------------------------------------------------------------------------
	// for (int i = 0; i < n_pts_; ++i) {
	// 	printf("[ ");
	// 	for (int j = 0; j < M_; ++j) {
	// 		printf("%llu ", hash_key_[i][j]);
	// 	}
	// 	printf("], id = %d\n", i + 1);
	// }
	// printf("\n");
}

// -----------------------------------------------------------------------------
bool M_SRP::calc_hash_code(			// calc hash code after random projection
	int   id,							// projection vector id
	const float *data)					// input data
{
	return calc_inner_product(dim_, proj_[id], data) >= 0 ? true : false;
}

// -----------------------------------------------------------------------------
uint64_t* M_SRP::compress_hash_code(// compress hash code with 64 bits
	bool *hash_code)					// input hash code
{
	uint64_t *hash_key = new uint64_t[M_];
	int shift = 0;
	for (int i = 0; i < M_; ++i) {
		uint64_t val = 0;
		int size = (i == M_-1 && m_%64 != 0) ? (m_ % 64) : 64;

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
int M_SRP::mcss(					// mcss
	int   top_k,						// top-k number of data objects
	int   check_k,						// check-k number of data objects
	const float *query,					// input query
	MaxK_List *list)					// top-k MC results (return)
{
	// -------------------------------------------------------------------------
	//  calculate the hash key (compressed hash code) of query
	// -------------------------------------------------------------------------
	bool *hash_code_q = new bool[m_];
	for (int i = 0; i < m_; ++i) {
		hash_code_q[i] = calc_hash_code(i, query);
	}
	uint64_t *hash_key_q = compress_hash_code(hash_code_q);

	// -------------------------------------------------------------------------
	//  find the candidates with largest matched values
	// -------------------------------------------------------------------------
	int total_bits = 64 * M_;
	int cand_size  = check_k + top_k;
	
	MaxK_List *check_list = new MaxK_List(cand_size);
	for (int i = 0; i < n_pts_; ++i) {
		uint32_t match = 0;
		for (int j = 0; j < M_; ++j) {
			match += table_lookup(hash_key_[i][j] ^ hash_key_q[j]);
			// match += num_bits(hash_key_[i][j] ^ hash_key_q[j]);
		}
		check_list->insert((float) (total_bits-match), i);
	}

	// -------------------------------------------------------------------------
	//  verification
	// -------------------------------------------------------------------------
	for (int i = 0; i < cand_size; ++i) {
		int id = check_list->ith_id(i);
		list->insert(calc_cosine(dim_, data_[id], query), id + 1);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] hash_code_q; hash_code_q = NULL;
	delete[] hash_key_q;  hash_key_q  = NULL;
	delete   check_list;  check_list  = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
uint32_t M_SRP::table_lookup(		// table lookup the match value
	uint64_t x)							// input uint64_t value
{
	return table16_[x & 0xffff] + table16_[(x>>16) & 0xffff] + 
		table16_[(x>>32) & 0xffff] + table16_[(x>>48) & 0xffff];
}

// -----------------------------------------------------------------------------
int M_SRP::num_bits(				// count number of 1 bits of input x
	uint64_t x)							// input uint64_t x value
{
	x = x - ((x >> 1) & 0x5555555555555555);
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
    x = ((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F);
    return (x*(0x0101010101010101)) >> 56;
}