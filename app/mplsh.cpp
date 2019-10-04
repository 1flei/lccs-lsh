#include "mplsh.h"

using namespace std;
using namespace MyCallbackRegister;

bool MPLSH_LSHKIT_REGISTED = registerCallback("mplsh_lshkit",
		"n qn d K L r dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = argAs<int>("n");
	int qn = argAs<int>("qn");
	int d = argAs<int>("d");
	int K = argAs<int>("K");
    int L = argAs<int>("L");
    // int T = argAs<int>("T");
    double r = argAs<double>("r");
	const float** data = argAs<const float**>("dataset");
	const float** query = argAs<const float**>("queryset");
	const Result** ground_truth = argAs<const Result**>("ground_truth");
	string output_filename = argAs<string>("output_filename");
    // int checked_candidate = argAs<int>("checked_candidate");


    // typedef MPLSH Index;
    typedef MPLSH_LSHKIT Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>(n, d, K, L, r);
        index->build(data);
        return index;
    };

    fprintf(fp.get(), "mplsh_lshkit  r=%f, K=%d, L=%d\n", r, K, L);
    // std::vector<int> nProbess = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    std::vector<int> nProbess = {1, 2, 4, 8, 16, 32, 64};

    const auto& fq = [&](Index &index, int k, int nProbes, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        index.query(nProbes*L, queryi, list);
    };

    benchmarkMinklist(qn, query, ground_truth, nProbess, fp.get(), fif, fq);
});