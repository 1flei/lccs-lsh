#ifndef __LCSB_H
#define __LCSB_H

#include <vector>
#include "../util.h"
#include <queue>
#include "../myTimer.h"
#include "../myndarray.h"

// -------------------------------------------------------------------------
struct Item {						// assis data structure for lcs matching
	// uint64_t *key_;						// hash code
	// int n_;								// total size
	int id_;							// data object id
	// int pos_;							// start bit position
};


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
		int m, 
        int step);

    // const std::vector<std::vector<uint64_t> > *datap;
    // void build(const std::vector<std::vector<uint64_t> > &data);
    void build(NDArray<2, uint64_t> &data);

	// -------------------------------------------------------------------------
	~LCSB();						// destructor

	// -------------------------------------------------------------------------
	int mcss(						// mcss
		int   check_k,					// check-k number of data objects
		const float *query,				// input query
		MaxK_List *list);				// top-k MC results  (return)

    template<typename F>
    // void forCandidates(int cand_size, const std::vector<uint64_t>& query, const F& f) 
    void forCandidates(int cand_size, const std::vector<uint64_t>& query, const F& f) 
    {  
        uint64_t *hash_key_q = (uint64_t *)&query[0];

        // -------------------------------------------------------------------------
        //  initialize the priority queue for candidate generation
        // -------------------------------------------------------------------------
        std::priority_queue<Elem> cand;
        int id    = -1;
        int loc   = -1;
        int left  = -1;
        int right = -1;
        int match = -1;

        // int step = (int) ceil(m_ * 0.05f);
        // int step = 1;
        MyTimer::pusht();
        for (loc = 0; loc < m_; loc = loc + step) {
            int pos = binary_search(loc, hash_key_q, match);

            if (match > 0) {
                // Elem obj(match, loc, pos-1, pos+1, list_[loc][pos].id_);
                Elem obj(match, loc, pos-1, pos+1, get_idx(list_[loc][pos]));
                cand.push(obj);
                
                // printf("i = %3d, pos = %5d, match = %d, idx=%d\n", loc, pos, match, get_idx(list_[loc][pos]));
            }
        }
        double t1 = MyTimer::popt();

        // -------------------------------------------------------------------------
        //  m-way merge sort and verification
        // -------------------------------------------------------------------------
        int count       = 0;
        int left_match  = -1;
        int right_match = -1;
        // int cand_size   = check_k + top_k;
        cand_size = cand_size > n_pts_ ? n_pts_ : cand_size;

        double tverify = 0;

        MyTimer::pusht();
        while (!cand.empty()) {
            Elem top = cand.top();
            id    = top.id_;
            loc   = top.loc_;
            left  = top.left_;
            right = top.right_;
            match = top.key_;

            cand.pop();

            printf("       top=%d, %d, %d, %d, %d\n", id, loc, left, right, match);

            // ---------------------------------------------------------------------
            //  verification
            // ---------------------------------------------------------------------
            if (checked_[id] < cnt_) {
                // MyTimer::pusht();
                // list->insert(calc_cosangle(dim_, data_[id], query), id + 1);
                f(id);
                checked_[id] = cnt_;
                // tverify += MyTimer::popt();
                if (++count >= cand_size) break;
            }

            // ---------------------------------------------------------------------
            //  get the next 2 possible candidates
            // ---------------------------------------------------------------------
            left_match = -1;			// case 1: left position
            if (left >= 0) cmp_hash_key(loc, left, hash_key_q, left_match);
            if (left_match == match) {
                // Elem obj(left_match, loc, left-1, right, list_[loc][left].id_);
                Elem obj(left_match, loc, left-1, right, get_idx(list_[loc][left]));
                cand.push(obj);
                continue;
            }

            right_match = -1;			// case 2: right position
            if (right < n_pts_) cmp_hash_key(loc, right, hash_key_q, right_match);
            if (right_match > left_match) {
                // Elem obj(right_match, loc, left, right+1, list_[loc][right].id_);
                Elem obj(right_match, loc, left, right+1, get_idx(list_[loc][right]));
                cand.push(obj);
            }
            else {
                // Elem obj(left_match, loc, left-1, right, list_[loc][left].id_);
                Elem obj(left_match, loc, left-1, right, get_idx(list_[loc][left]));
                cand.push(obj);
            }
        }
        double t2 = MyTimer::popt();
        printf("    t1, t2 =%f, %f\n", t1, t2);
	    // MyTimer::printAll();
        // printf("match_cnt=%d\n", MATCH_CNT);
        // printf("\n");

        // -------------------------------------------------------------------------
        //  update <cnt_> and <checked_>
        // -------------------------------------------------------------------------
        cnt_++;
        if (cnt_ == 4294967295U) {
            cnt_ = 1;
            memset(checked_, 0, sizeof(uint32_t) * n_pts_);
        }
    }

protected:



	//O(dim * 2^K + dim * n), which is accectable if K is chosen carefully 
	void init_extent();
	// static const int K = 16;
	// static const int TableSize = (1 << K);
	// std::vector<std::array<int32_t, TableSize + 1> > extents;
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

    inline int get_idx(const Item& x) const
    {
        return x.id_;
    }
    inline uint64_t* get_key(const Item& x) const
    {
        return hash_key_[x.id_];
        // return x.key_;
    }

	int   n_pts_;					// cardinality of dataset
	int   m_;						// number of hash functions
    int   logn;
    int   step;

	float    **proj_;				// random projection vectors
	Item     **list_;				// sorted list of cmpr hash key of all objs 
	uint64_t **hash_key_;			// compressed hash key of all objs
	uint8_t  *table16_;				// lookup table for uint16_t value

	uint32_t cnt_;					// counter for checked ID
	uint32_t *checked_;				// checked ID for each query


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


    // -----------------------------------------------------------------------------
    int cmp_item(						// item compare func for qsort (ascending)
        const Item& e1,						// 1st element
        const Item& e2, 
        int pos); 					// 2nd element

	// -------------------------------------------------------------------------
	inline int table_lookup(		// table lookup the match value
		uint64_t x);					// input uint64_t value
};

#endif // __LCSB_H
