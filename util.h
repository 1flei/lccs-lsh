#pragma once

#include <memory>
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstdint> 
#include <vector>
#include <random>

#include <unistd.h>
#include "def.h"
#include "pri_queue.h"


// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
//  uitlity functions
// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_data(						// read data set from disk
	int   n,							// number of data points
	int   d,							// dimensionality
	const char *fname,					// address of data
	Scalar **data);						// data (return)

// -----------------------------------------------------------------------------
int read_data_binary(						// read data set from disk
	int   n,							// number of data points
	int   d,							// dimensionality
	const char *fname,					// address of data
	Scalar **data);						// data (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int    qn,							// number of query objects
	const  char *fname,					// address of truth set
	Result **R);		
	
// -----------------------------------------------------------------------------
Scalar calc_inner_product(			// calc inner product
	int   dim,							// dimension
	const Scalar *p1,					// 1st point
	const Scalar *p2);					// 2nd point
	

// -----------------------------------------------------------------------------
Scalar calc_l1_dist(					// calc L1 distance
	int   dim,							// dimension
	const Scalar *p1,					// 1st point
	const Scalar *p2);					// 2nd point

// -----------------------------------------------------------------------------
Scalar calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
int get_hits(						// get the number of hits between two ID list
	int   k,							// top-k value
	int   t,							// top-t value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
Scalar calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results
	MinK_List *list);					// results returned by algorithms
// -----------------------------------------------------------------------------
Scalar calc_weighted_dist2(			// calc inner product
	int   dim,							// dimension
	const Scalar *w,
	const Scalar *p1,					// 1st point
	const Scalar *p2);					// 2nd point

inline int get_num_bits8(uint8_t x)				//get the number of 1 in the binary representation of u
{
	x = (x&0x55) + ((x>>1)&0x55);
	x = (x&0x33) + ((x>>2)&0x33);
	x = (x&0x0f) + ((x>>4)&0x0f);
	return x;
}

inline int get_num_bits64(uint64_t x)			////get the number of 1 in the binary representation of x
{
	x = x - ((x >> 1) & 0x5555555555555555);
    x = (x & 0x3333333333333333) +
        ((x >> 2) & 0x3333333333333333);
    x = ((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F);
    return (x*(0x0101010101010101))>>56;
}

const int PrefixTableSize = 1 << 16;
extern std::array<uint8_t, PrefixTableSize> _prefix_table;

inline bool init_prefix_table()
{
	for (int i = 0; i < PrefixTableSize; i++) {
		//calculate the prefix-1 of i, since it will be run only once, implement using stupid way
		_prefix_table[i] = 0;
		for (int j = 0; j < 16; j++) {
			int mask = 1 << (15 - j);
			if (i&mask) {
				_prefix_table[i]++;
			}
			else {
				break;
			}
		}
	}
	return true;
}

inline int get_num_prefix(uint16_t u)
{
	static bool initialized = init_prefix_table();
	return _prefix_table[u];
}

inline int get_num_prefix(uint32_t u)
{
	int a = get_num_prefix(uint16_t(u >> 16));
	if (a != 16) {
		return a;
	}
	int b = get_num_prefix(uint16_t(u & 0xffff));
	return a + b;
}

inline int get_num_prefix(uint64_t u)
{
	int a = get_num_prefix(uint16_t(u >> 48));
	if (a != 16) {
		return a;
	}
	int b = get_num_prefix(uint16_t((u >> 32) & 0xffff));
	if (b != 16) {
		return a + b;
	}
	int c = get_num_prefix(uint16_t((u >> 16) & 0xffff));
	if (c != 16) {
		return a + b + c;
	}
	int d = get_num_prefix(uint16_t(u & 0xffff));
	return a + b + c + d;
}



inline uint64_t hash_combine(uint64_t h0, uint64_t h1)
{
    return h0 ^ (0x9e3779b9 + (h0<<6) + (h0>>2) + h1);
}

inline int log2i(unsigned x) {
    return sizeof(int) * 8 - __builtin_clz(x) - 1;
}
inline int log2ll(unsigned long long x) {
    return sizeof(int) * 8 - __builtin_clzll(x) - 1;
}



int calc_hamming_dist(			// calc inner product
	int   dim,		
	const uint8_t *p1,					// 1st point
	const uint8_t *p2);					// 2nd point

int calc_hamming_dist(			// calc inner product
	int   dim,		
	const uint64_t *p1,					// 1st point
	const uint64_t *p2);					// 2nd point

Scalar calc_angle(				// calc angle
	int   dim,							// dimension
	const Scalar *p1,					// 1st point
	const Scalar *p2);					// 2nd point


Scalar calc_cosangle(				// calc cos(angle)
	int   dim,							// dimension
	const Scalar *p1,					// 1st point
	const Scalar *p2);					// 2nd point

Scalar calc_ratio(
	int k, 
	const Result *Rs, 
	MinK_List *list);

Scalar calc_ratio(
	int k, 
	const Result *Rs, 
	MaxK_List *list);

#if __cplusplus < 201402L
//make_unique, which is available in c++17
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#else
using std::make_unique;
#endif



// -----------------------------------------------------------------------------
template<class ScalarType>
ScalarType calc_l2_sqr(					// calc L2 square distance
	int   dim,							// dimension
	const ScalarType *p1,					// 1st point
	const ScalarType *p2)					// 2nd point
{
	ScalarType diff(0);
	ScalarType ret(0);
	for (int i = 0; i < dim; ++i) {
		diff = p1[i] - p2[i];
		ret += diff * diff;
	}
	return ret;
}

// -----------------------------------------------------------------------------
template<class ScalarType>
ScalarType calc_l2_dist(					// calc L2 distance
	int   dim,							// dimension
	const ScalarType *p1,					// 1st point
	const ScalarType *p2)					// 2nd point
{
	return sqrt(calc_l2_sqr(dim, p1, p2));
}


template<class Iter>
void printVec(const Iter &its, int dim)
{
	Iter it = its;
	for(int i=0;i<dim;i++){
		std::cout << *it;
		++it;
		if(i==dim-1){
			std::cout << std::endl;
		} else{
			std::cout << ", ";
		}
	}
}

template<class Iter>
void printVec(const Iter& its, const Iter& ite, const std::string &header=""){
//	if(verbosity){
	std::cout << header;
	for(Iter it=its; ; ){
		std::cout << *it;
		++it;
		if(it!=ite){
			std::cout << ", ";
		} else{
			std::cout << std::endl;
			break;
		}
	}
//	}
}
template<class Iter, class OStream>
void printVec(const Iter& its, const Iter& ite, const std::string &header="", OStream &os=std::cout){
//	if(verbosity){
	os << header;
	for(Iter it=its; ; ){
		os << *it;
		++it;
		if(it!=ite){
			os << ", ";
		} else{
			os << std::endl;
			break;
		}
	}
//	}
}

//return idx such that xs[idx[i]] is sorted for i in range(beign, end)
//require begin and end are randomly accessable!!
template<typename T, typename F>
std::vector<int> argsort(const T& begin, const T& end, const F& cmp)
{
	size_t len = distance(begin, end);
	std::vector<int> idx(len);
	for(int i=0;i<idx.size();i++){
		idx[i] = i;
	}
	std::sort(idx.begin(), idx.end(), [&](int a, int b){
		return cmp(*(begin+a), *(begin+b));
	});
	return idx;
}
template<typename T>
std::vector<int> argsort(const T& begin, const T& end)
{
	size_t len = distance(begin, end);
	std::vector<int> idx(len);
	for(int i=0;i<idx.size();i++){
		idx[i] = i;
	}
	std::sort(idx.begin(), idx.end(), [&](int a, int b){
		return *(begin+a)<*(begin+b);
	});
	return idx;
}

template<typename Iter>
void scatter(const Iter& begin, std::vector<int> &idx)
{
    using T = typename std::iterator_traits<Iter>::value_type;
	std::vector<T> tmpxs(begin, begin+idx.size());
	for(int i=0;i<idx.size();i++){
		*(begin+i) = tmpxs[idx[i]];
	}
}

template<typename T> 
std::vector<int> getExtent(std::vector<T>& toSplit, std::vector<T>& refs)
{
	std::vector<int> ret;
	ret.reserve(refs.size()+1);

	int curLoc = 0;
	for(T& r:refs){
		while(curLoc<toSplit.size() && r>toSplit[curLoc]){
			curLoc++;
		}
		ret.push_back(curLoc);
	}
	ret.push_back(toSplit.size());
	return ret;
}

template<class uintt>
struct CountMarkerU
{
	std::vector<uintt> markCount;
	uintt curCnt;
	CountMarkerU(int sz=0):markCount(sz), curCnt(1){}


	void resize(int n){
		markCount.resize(n);
		fill(markCount.begin(), markCount.end(), 0);
	}

	void mark(int n){
		markCount[n] = curCnt;
	}
	bool isMarked(int n){
		return markCount[n] >= curCnt;
	}
	void clear(){
		if(curCnt==~uintt(0)){
			curCnt=1;
			markCount.clear();
		} else{
			curCnt++;
		}
	}
};

using CountMarker = CountMarkerU<unsigned>;

//array version of CountMatker
template<class T> 
struct CountArray
{
	CountMarker marker;
	CountArray(int sz):count(sz), marker(sz){}

	std::vector<T> count;

	void resize(int n){
		marker.resize(n);
		count.resize(n);
		marker.clear();
	}

	T& operator[](int i){
		if(!marker.isMarked(i)){
			marker.mark(i);
			count[i] = 0;
			return count[i];
		}
		return count[i];
	}

	void clear(){
		marker.clear();
	}
};


//combine hash signature using pertubation hash
template<class INTTYPE>
struct HashCombinator
{
	HashCombinator(int n, uint64_t range=0xffffffffffffffff)
		: n(n)
	{
		std::mt19937 rng(time(NULL));
		if(range < 0x7fffffffffffffff){
			mask = 1;
			while(mask < range){
				mask <<= 1;
			}
			--mask;
		} else{
			mask = 0xffffffffffffffff;
		}

		keys.resize(n*SIZEPERONE);
		for(int i=0;i<n*SIZEPERONE;i++){
			for(int j=0;j<256;j++){
				uint64_t r = rng();
				keys[i][j] = r & mask;
			}
		}
	}

	static const int SIZEPERONE = sizeof(INTTYPE);
	int n;
	uint64_t mask;

	std::vector<std::array<uint64_t, 256> > keys;

	template<typename IT> 
	uint64_t hash_combine(const IT& beg)
	{
		uint64_t ret = 0;
		IT cur = beg;
		for(int i=0;i<n;i++){
			INTTYPE cur_hash = *cur;
			for(int j=0;j<SIZEPERONE;j++){
				uint8_t kij = cur_hash & 0xff;
				ret ^= keys[i*SIZEPERONE+j][kij];
				cur_hash >>= 8;
			}
			++cur;
		}
		return ret;
	}
};