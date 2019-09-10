#include "mplsh.h"

using namespace std;
using namespace MyCallbackRegister;

bool MYMPLSH_REGISTED = registerCallback("mymplsh",
		"n qn d K L r dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int L = algAs<int>("L");
    // int T = algAs<int>("T");
    double r = algAs<double>("r");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef MPLSH Index;
    typedef MPLSH Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>(n, d, K, L, r);
        index->build(data);        
        return index;
    };

    fprintf(fp.get(), "mymplsh  r=%f, K=%d, L=%d\n", r, K, L);
    std::vector<int> nProbess = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

    const auto& fq = [&](Index &index, int k, int nProbes, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        index.query(nProbes, queryi, list);
    };

    benchmarkMinklist(qn, query, ground_truth, nProbess, fp.get(), fif, fq);
});

bool MPLSH_LSHKIT_REGISTED = registerCallback("mplsh_lshkit",
		"n qn d K L r dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int L = algAs<int>("L");
    // int T = algAs<int>("T");
    double r = algAs<double>("r");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef MPLSH Index;
    typedef MPLSH_LSHKIT Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>(n, d, K, L, r);
        index->build(data);
        return index;
    };

    fprintf(fp.get(), "mplsh_lshkit  r=%f, K=%d, L=%d\n", r, K, L);
    std::vector<int> nProbess = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

    const auto& fq = [&](Index &index, int k, int nProbes, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        index.query(nProbes, queryi, list);
    };

    benchmarkMinklist(qn, query, ground_truth, nProbess, fp.get(), fif, fq);
});