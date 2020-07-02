#pragma once
#include "def.h"

// -----------------------------------------------------------------------------
//  struct Result
// -----------------------------------------------------------------------------
struct Result {						// structure for furthest neighbor / hash value
	Scalar key_;							// distance / random projection value
	int   id_;							// object id

	Result() : key_(0.), id_(-1) {}
	Result(Scalar k, int id) : key_(k), id_(id) {}

};