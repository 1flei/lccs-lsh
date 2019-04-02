#include "rh_lcs_sort.h"

using namespace std;

bool RH_LCS_SORT::registered = MyCallbackRegister::registerCallback("rh_lcs_sort", 
		"n qn d K M dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int M = algAs<int>("M");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");

    const auto fif = [&](){
		auto lsh = make_unique<RH_LCS_SORT>(n, d, K, data, M);
        return lsh;
    };

	const auto fq = [&](RH_LCS_SORT &lsh, int k, const float* queryi, MinK_List* list) {
        lsh.kmc_angle(k, queryi, list);
    };

    benchmarkMinklist(qn, d, query, ground_truth, output_filename, fif, fq);
});

// -----------------------------------------------------------------------------
RH_LCS_SORT::RH_LCS_SORT(					// constructor
	int n,								// cardinality of dataset
	int d,								// dimensionality of dataset
	int K,								// number of hash tables
	const float **data, 
    int M)					// data objects
    : M_(M)
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	K_     = K;
	data_  = data;

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
    pivots.resize(K);
    for(int i=0;i<K;i++){
        pivots[i] = rand()%n;
    }

    for(int i=0;i<K;i++){
        for(int j=i+1;j<K;j++){
            orders.emplace_back(i, j);
        }
    }
    random_shuffle(orders.begin(), orders.end());

	bulkload();
}

// -----------------------------------------------------------------------------
RH_LCS_SORT::~RH_LCS_SORT()					// destructor
{
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void RH_LCS_SORT::bulkload()			// bulkloading
{
    codes.reserve(n_pts_);
	for (int i = 0; i < n_pts_; ++i) {
        codes.emplace_back(get_proj_vector(data_[i]));
	}
    samIndexer = make_unique<SALCS>(codes[0].size(), codes);
//	printf("finish computing hash-sigs\n");
}

std::vector<int> RH_LCS_SORT::get_proj_vector(const float *data)
{
    std::vector<int> ret;
    std::vector<Scalar> dists(K_);
	for (int i = 0; i < K_; ++i) {
        dists[i] = calc_cosangle(dim_, data, data_[pivots[i]]);
    }

    int r = 0;
    int m = 0;
    int sig = 0;
    for(auto& uv:orders){
        int u = uv.first;
        int v = uv.second;

        if(dists[u] < dists[v]){
            sig = (sig<<1) +1;
        } else{
            sig = (sig<<1);
        }
        if(++m==M_){
            ret.push_back((u*K_+v+1)*(1<<M_) + sig);
            m = 0;
            sig = 0;
        }
    }
    return ret;
}


// -----------------------------------------------------------------------------
int RH_LCS_SORT::kmc_angle(					// c-k-AMC search
	int   top_k,						// top-k value
	const float *query,					// input query
	MinK_List *list)					// top-k MC results (return)
{
//	bool *mc_query = new bool[K_];
//	get_proj_vector(query, mc_query);
	auto codes = get_proj_vector(query);
    int nCandidates = top_k+K_+K_*(K_-1)/2/M_/dim_*10;

    const auto& f = [&](int idx){
        float angle = calc_angle(dim_, data_[idx], query);
		list->insert(angle, idx + 1);
    };
    samIndexer->klcs(codes, nCandidates, f);

	return 0;
}
