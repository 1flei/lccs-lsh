#ifndef __M_SRP_H
#define __M_SRP_H

// -----------------------------------------------------------------------------
//  m-Sign-Random Projection LSH (or simply m-srp) is used to solve the 
//  problem of Maximum Cosine Similarity Search (MCSS)
// -----------------------------------------------------------------------------
class M_SRP {
public:
	M_SRP(							// constructor
		int   n,						// cardinality of dataset
		int   d,						// dimensionality of dataset
		int   m,						// number of hash functions
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~M_SRP();						// destructor

	// -------------------------------------------------------------------------
	int mcss(						// mcss
		int   top_k,					// top-k number of data objects
		int   check_k,					// check-k number of data objects
		const float *query,				// input query
		MaxK_List *list);				// top-k MC results  (return)

protected:
	int   n_pts_;					// cardinality of dataset
	int   dim_;						// dimensionality of dataset
	int   m_;						// number of hash functions
	int   M_;						// number of compress uint64_t hash code
	const float **data_;			// data objects

	float    **proj_;				// random projection vectors
	uint64_t **hash_key_;			// hash key (compressed hash code)
	uint32_t *table16_;				// table to record the number of 1 bits 

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
		bool *hash_code);				// input hash code
	
	// -------------------------------------------------------------------------
	uint32_t table_lookup(			// table lookup the match value
		uint64_t x);					// input uint64_t value
	
	// -------------------------------------------------------------------------
	int num_bits(					// count number of 1 bits of input x
		uint64_t x);					// input uint64_t x value
};

#endif // __M_SRP_H
