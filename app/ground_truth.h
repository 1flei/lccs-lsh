#pragma once

#include <memory>
#include "../myTimer.h"
#include "../util.h"
#include "../pri_queue.h"
#include "../def.h"
#include "../register.h"

extern bool ground_truth_registed;
extern bool ground_truth_for_ws_registed;
extern bool ground_truth_furthest_registed;
extern bool ground_truth_angle_registed;

// -----------------------------------------------------------------------------
//  interface of this package
// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set);			// address of truth set

int ground_truth_for_weighted_space(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const float **weight,				// query set
	const char  *truth_set);				// address of truth set

int ground_truth_furthest(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set);				// address of truth set

int ground_truth_angle(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set);				// address of truth set