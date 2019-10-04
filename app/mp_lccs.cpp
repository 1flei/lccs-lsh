#include "mp_lccs.h"

using namespace std;
using namespace MyCallbackRegister;

bool MP_LCCS_REGISTED = registerCallback("mp_lccs",
		"n qn d L r p dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
    int L = algAs<int>("L");
    double r = algAs<double>("r");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");

    int p = algAs<int>("p");
    int defaultLookBack=7;

    // typedef MPLSH Index;
    typedef mylccs::MP_LCCS_L2 Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>(n, d, L, r, defaultLookBack, p);
        index->build(data);
        return index;
    };

    fprintf(fp.get(), "mplsh_lshkit  r=%f, L=%d\n", r, L);
    // std::vector<int> nProbess = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    std::vector<int> nScanss = {1, 2, 4, 8, 16, 32};

    const auto fq = [&](Index &index, int k, int nScans, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        const auto f = [&](int idx){
            float l2dist = calc_l2_dist(d, data[idx], queryi);
            list->insert(l2dist, idx + 1);
        };
        index.query(nScans, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, nScanss, fp.get(), fif, fq);
});