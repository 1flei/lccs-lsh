#include "srp_lcs.h"

using namespace std;

// bool SRP_LCS::registered = MyCallbackRegister::registerCallback("srp_lcs", 
// 		"n qn d K dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
// 	using namespace MyCallbackRegister;
// 	int n = algAs<int>("n");
// 	int qn = algAs<int>("qn");
// 	int d = algAs<int>("d");
// 	int K = algAs<int>("K");
// 	const float** data = algAs<const float**>("dataset");
// 	const float** query = algAs<const float**>("queryset");
// 	const Result** ground_truth = algAs<const Result**>("ground_truth");
// 	string output_filename = algAs<string>("output_filename");

//     const auto fif = [&](){
// 		auto lsh = make_unique<SRP_LCS>(n, d, K, data);
//         return lsh;
//     };

// 	const auto fq = [&](SRP_LCS &lsh, int k, const float* queryi, MinK_List* list) {
//         lsh.kmc_angle(k, queryi, list);
//     };

//     benchmarkMinklist(qn, query, ground_truth, output_filename, fif, fq);
// });

// -----------------------------------------------------------------------------
SRP_LCS::SRP_LCS(					// constructor
	int n,								// cardinality of dataset
	int d,								// dimensionality of dataset
	int K,								// number of hash tables
	const float **data, 
    int M)					// data objects
    : samIndexer(K), M_(M)
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
	gen_random_vectors();
	bulkload();
}

// -----------------------------------------------------------------------------
SRP_LCS::~SRP_LCS()					// destructor
{
	if (proj_ != NULL) {
		for (int i = 0; i < K_; ++i) {
			delete[] proj_[i];	proj_[i] = NULL;
		}
		delete[] proj_;	proj_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void SRP_LCS::gen_random_vectors()	// generate random projection vectors
{
    std::normal_distribution<double> gaussian(0.);
    std::random_device rd;
    std::default_random_engine rng(rd() );

	proj_ = new float*[K_*M_];
	for (int i = 0; i < K_; ++i) {
		proj_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			proj_[i][j] = gaussian(rng);
		}
	}
}

// -----------------------------------------------------------------------------
void SRP_LCS::bulkload()			// bulkloading
{
	for (int i = 0; i < n_pts_; ++i) {
        auto codes = get_proj_vector(data_[i]);
        samIndexer.insert(codes);
	}
//	printf("finish computing hash-sigs\n");
}


std::vector<int> SRP_LCS::get_proj_vector(const float *data)
{
    std::vector<int> ret;
	for (int i = 0; i < K_; ++i) {
        int sigi = 0;
        for(int j=0;j<M_;j++){
    		float sum = calc_inner_product(dim_, proj_[i], data);
            if(sum >= 0){
                sigi = (sigi<<1) +1;
            } else{
                sigi = sigi<<1;
            }
        }
        ret.push_back( (i+1)*(1<<M_) + sigi);
	}
    return ret;
}


// -----------------------------------------------------------------------------
int SRP_LCS::kmc_angle(					// c-k-AMC search
	int   top_k,						// top-k value
	const float *query,					// input query
	MinK_List *list)					// top-k MC results (return)
{
//	bool *mc_query = new bool[K_];
//	get_proj_vector(query, mc_query);
	auto codes = get_proj_vector(query);
    int nCandidates = top_k+K_*M_;

    const auto& f = [&](int idx){
        float angle = calc_angle(dim_, data_[idx], query);
		list->insert(angle, idx + 1);
    };
    samIndexer.klcs(codes, nCandidates, f);

	return 0;
}
