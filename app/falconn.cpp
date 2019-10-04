#include "falconn.h"

using namespace std;
using namespace MyCallbackRegister;

bool FALCONN_REGISTED = registerCallback("falconn",
		"n qn d L cp_dim dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	// int nBits = algAs<int>("nHashBits");
    int L = algAs<int>("L");
    int K = algAs<int>("K");

    int lastCpDim = algAs<int>("cp_dim");

	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef MPLSH Index;
    typedef FALCONN_INDEX Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>(n, d, L, K, lastCpDim);
        index->build(data);
        return index;
    };

    fprintf(fp.get(), "falconn  K=%d, L=%d\n", K, L);
    std::vector<int> nProbess = {1, 2, 4, 8, 16, 32};

    const auto& fq = [&](Index &index, int k, int nProbes, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        index.query(nProbes*L, queryi, list);
    };

    benchmarkMinklist(qn, query, ground_truth, nProbess, fp.get(), fif, fq);
});