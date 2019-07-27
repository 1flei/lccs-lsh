#ifndef __INTERFACE_H
#define __INTERFACE_H

// -----------------------------------------------------------------------------
//  interface of maximum cosine similarity search (MCSS, or minimum angle)
// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth MCSS results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set);			// address of truth set

// -----------------------------------------------------------------------------
int linear_scan(					// find top-k mip using linear_scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder); 		// output folder

// -----------------------------------------------------------------------------
int m_srp(							// m-srp for MCSS
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// number of projections (concatenation)
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder);		// output folder

// -----------------------------------------------------------------------------
int kl_srp(							// (k,l)-srp for MCSS
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   k,							// number of projections (concatenation)
	int   l,							// number of hash tables (repetion)
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder);		// output folder

// -----------------------------------------------------------------------------
int lcsb(							// lcsb for MCSS
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// number of projections (concatenation)
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder);		// output folder

#endif // __INTERFACE_H
