#ifndef __KL_SRP_H
#define __KL_SRP_H

#include <vector>
#include <unordered_map>

using namespace std;

// -----------------------------------------------------------------------------
//  (k,l)-Sign-Random Projection LSH (or simply (k,l)-srp) is used to solve the 
//  problem of Maximum Cosine Similarity Search (MCSS)
// -----------------------------------------------------------------------------
class KL_SRP {
public:
	KL_SRP(							// constructor
		int n,							// cardinality of dataset
		int d,							// dimensionality of dataset
		int k,							// number of projections per table
		int l,							// number of hash tables
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~KL_SRP();						// destructor

	// -------------------------------------------------------------------------
	int mcss(						// mcss
		int   top_k,					// top-k number of data objects
		int   check_k,					// check-k number of data objects
		const float *query,				// input query
		MaxK_List *list);				// top-k MC results  (return)

protected:
	int   n_pts_;					// cardinality of dataset
	int   dim_;						// dimensionality of dataset
	int   k_;						// number of projections per table
	int   l_;						// number of hash tables
	int   m_;						// number of standard hash functions
	const float **data_;			// data objects

	vector<unordered_map<int, vector<int> > > tables_; // hash tables
	float ***proj_;					// random projection vectors, O(l*k*d)
	uint32_t *standard_hash_;		// standard hash function

	uint32_t cnt_;					// counter for checked ID
	uint32_t *checked_;				// checked ID for each query
	
	// -------------------------------------------------------------------------
	void gen_random_vectors();		// generate random projection vectors

	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading

	// -------------------------------------------------------------------------
	int calc_hash_value(			// calculate hash value of input data
        int   hash_id,                  // hash table id
		const float *data);				// input data
};

#endif // __KL_SRP_H
