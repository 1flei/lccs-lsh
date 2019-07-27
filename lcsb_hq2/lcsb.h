#ifndef __LCSB_H
#define __LCSB_H

#include <vector>
#include "../myTimer.h"
using namespace std;

// -------------------------------------------------------------------------
struct Item {						// assis data structure for lcs matching
	uint64_t *key_;						// hash code
	int n_;								// total size
	int id_;							// data object id
	int pos_;							// start bit position
};

// -------------------------------------------------------------------------
int cmp_item(						// item cmp func for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
struct Elem {						// assistant data structure for mcss
    int key_;							// match size of the hash code
	int loc_;							// location of the hash code
	int left_;							// left position of the bucket list
	int right_;							// right position of the bucket list
    int id_;							// data object id

	Elem() {
		key_   = -1;
		loc_   = -1;
		left_  = -1;
		right_ = -1;
		id_    = -1;
	}

	Elem(int key, int loc, int left, int right, int id) {
		key_   = key;
		loc_   = loc;
		left_  = left;
		right_ = right;
		id_    = id;
	}

	bool operator<(const Elem &e) const { 
		return key_ < e.key_; 
	}
};

// -----------------------------------------------------------------------------
//  Longest Common Substring Bucketing (or simply LCSB) is used to solve the 
//  problem of Maximum Cosine Similarity Search (MCSS)
// -----------------------------------------------------------------------------
class LCSB {
public:
	LCSB(							// constructor
		int n,							// cardinality of dataset
		int d,							// dimensionality of dataset
		int m,							// number of hash functions
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~LCSB();						// destructor

	// -------------------------------------------------------------------------
	int mcss(						// mcss
		int   top_k,					// top-k number of data objects
		int   check_k,					// check-k number of data objects
		const float *query,				// input query
		MaxK_List *list);				// top-k MC results  (return)



	//O(dim * 2^K + dim * n), which is accectable if K is chosen carefully 
	void init_extent();
	// static const int K = 12;
	// static const int TableSize = (1 << K);
	std::vector<std::vector<int32_t> > extents;

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
		int loc64p = (loc64 + 1) % (m_/64);
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

protected:
	int   n_pts_;					// cardinality of dataset
	int   dim_;						// dimensionality of dataset
	int   m_;						// number of hash functions
	const float **data_;			// data objects
	int  logn;

	float    **proj_;				// random projection vectors
	Item     **list_;				// sorted list of cmpr hash key of all objs 
	uint64_t **hash_key_;			// compressed hash key of all objs
	uint8_t  *table16_;				// lookup table for uint16_t value

	uint32_t cnt_;					// counter for checked ID
	uint32_t *checked_;				// checked ID for each query

	// -------------------------------------------------------------------------
	void gen_random_vectors();		// generate random projection vectors

	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading

	// -------------------------------------------------------------------------
	bool calc_hash_code(			// calc hash code after random projection
		int   id,						// projection vector id
		const float *data);				// input data

	// -------------------------------------------------------------------------
	uint64_t* compress_hash_code(	// compress hash code with 64 bits
		int  n,							// compress length
		bool *hash_code);				// input hash code

	// -------------------------------------------------------------------------
	int binary_search(				// binary search the loc of input hash key
		int loc,						// location of the hash code
		uint64_t *hash_key,				// input hash key
		int &match);					// number of match item (return)
		
	// -------------------------------------------------------------------------
	int cmp_hash_key(				// cmp func for input hash key
		int loc,						// location of hash key
		int pos,						// position of hash key
		uint64_t *hash_key,				// input hash key
		int &match);					// number of match bits (return)

	// -------------------------------------------------------------------------
	inline int table_lookup(		// table lookup the match value
		uint64_t x);					// input uint64_t value
};

#endif // __LCSB_H
