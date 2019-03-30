#pragma once
#include "def.h"

// -----------------------------------------------------------------------------
//  struct Result
// -----------------------------------------------------------------------------
struct Result {						// structure for furthest neighbor / hash value
	Scalar key_;							// distance / random projection value
	int   id_;							// object id
};