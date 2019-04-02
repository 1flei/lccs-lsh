#pragma once 

#include "boost/dynamic_bitset.hpp"
#include <vector>
#include "../util.h"
#include "../pri_queue.h"
#include <random>
#include "sort_lcs.h"
#include "../register.h"
#include "../benchmark_util.h"

// -----------------------------------------------------------------------------
//  Sign-Random Projection LSH (RH_LCS_SORT) is used to solve the problem of 
//  c-Approximate Maximum Cosine (c-AMC) search
// 
//  the idea was introduced by Moses S Charikar in his paper "Similarity 
//  estimation techniques from rounding algorithms", In Proceedings of the 
//  thiry-fourth annual ACM symposium on Theory of computing (STOC), pages 
//  380–388, 2002.
// -----------------------------------------------------------------------------
class RH_LCS_SORT {
public:
	RH_LCS_SORT(
		int n,							// cardinality of dataset
		int d,							// dimensionality of dataset
		int K,							// number of hash tables
		const float **data, 
		int M = 4);			// data objects

	// -------------------------------------------------------------------------
	~RH_LCS_SORT();						// destructor

	// -------------------------------------------------------------------------
	int kmc(						// c-k-AMC search
		int   top_k,					// top-k value
		const float *query,				// input query
		MaxK_List *list);				// top-k MC results  (return)

	int kmc_angle(						// c-k-AMC search
		int   top_k,					// top-k value
		const float *query,				// input query
		MinK_List *list);				// top-k MC results  (return)

    static bool registered;
protected:
	int   n_pts_;					// cardinality of dataset
	int   dim_;						// dimensionality of dataset
	int   K_;						// number of hash tables
	const float **data_;			// data objects
	int   M_; 						// the number of hash functions per signature

//	bool  **hash_code_;				// hash code of data objects

	// float **proj_;					// random projection vectors
    std::vector<int> pivots;
    std::vector<std::pair<int, int> > orders;

	// -------------------------------------------------------------------------
	void gen_random_vectors();		// generate random projection vectors

	// -------------------------------------------------------------------------
	void bulkload();				// bulkloading

	std::vector<int> get_proj_vector(const float *data);

    std::vector<std::vector<int> > codes;

	std::unique_ptr<SALCS> samIndexer;
};