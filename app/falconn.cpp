#include "falconn.h"

using namespace std;
using namespace MyCallbackRegister;

bool FALCONN_REGISTED = registerCallback("falconn",
		"n qn d L dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = argAs<int>("n");
	int qn = argAs<int>("qn");
	int d = argAs<int>("d");
	// int nBits = argAs<int>("nHashBits");
    int L = argAs<int>("L");
    int K = argAs<int>("K");

	const float** data = argAs<const float**>("dataset");
	const float** query = argAs<const float**>("queryset");
	const Result** ground_truth = argAs<const Result**>("ground_truth");
	string output_filename = argAs<string>("output_filename");
    // int checked_candidate = argAs<int>("checked_candidate");

    bool useBitPackedHashTable = hasArg("bit_packed");

    int numRotation = 1;

    // typedef MPLSH Index;
    typedef FALCONN_INDEX Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>(n, d, L, K, numRotation, useBitPackedHashTable);
        index->build(data);
        return index;
    };

    fprintf(fp.get(), "falconn  K=%d, L=%d\n", K, L);
    std::vector<int> nProbess = {0, 1, 2, 4, 8, 16, 32};

    const auto& fq = [&](Index &index, int k, int nProbes, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        index.query(nProbes*L, queryi, list);
    };

    benchmarkMinklist(qn, query, ground_truth, nProbess, fp.get(), fif, fq);
});