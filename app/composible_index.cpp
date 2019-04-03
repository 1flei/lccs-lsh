#include "composible_index.h"

#include "../hashAlg/srp.h"
#include "../hashAlg/pivots.h"
#include "../bucketAlg/hamming.h"
#include "../bucketAlg/lcs.h"

using namespace std;
using namespace MyCallbackRegister;

bool SRPP_LINEAR_SCAN_REGISTED = registerCallback("srppair_scan", 
		"n qn d K dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int M = sizeof(SigType)*8;
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");


    typedef ComposibleIndex<SRPPair, HammingLinearScan> SRPP_HLS;

    const auto& fif = [&](){
		auto index = make_unique<SRPP_HLS>();
        index->initHasher(d, K, M);
        index->initBucketer(K);
        index->build(n, data);
        return index;
    };

	const auto& fq = [&](SRPP_HLS &index, int k, const float* queryi, MinK_List* list) {
        int nCandidates = k;
        // int nCandidates = n;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, d, query, ground_truth, output_filename, fif, fq);
});

bool SRP_LINEAR_SCAN_REGISTED = registerCallback("srp_scan", 
		"n qn d K dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int M = sizeof(SigType)*8;
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");


    typedef ComposibleIndex<SRP, HammingLinearScan> SRP_HLS;

    const auto& fif = [&](){
		auto index = make_unique<SRP_HLS>();
        index->initHasher(d, K, M);
        index->initBucketer(K);
        index->build(n, data);
        return index;
    };

	const auto& fq = [&](SRP_HLS &index, int k, const float* queryi, MinK_List* list) {
        int nCandidates = k*10;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, d, query, ground_truth, output_filename, fif, fq);
});

bool PIVOTSBINARY_LINEAR_SCAN_REGISTED = registerCallback("pb_scan", 
		"n qn d K dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int M = sizeof(SigType)*8;
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");


    typedef ComposibleIndex<PivotsBinary, HammingLinearScan> PB_HLS;

    const auto& fif = [&](){
		auto index = make_unique<PB_HLS>();
        index->initHasher(d, K, M, n, data);
        index->initBucketer(K);
        index->build(n, data);
        return index;
    };

	const auto& fq = [&](PB_HLS &index, int k, const float* queryi, MinK_List* list) {
        int nCandidates = k*10;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, d, query, ground_truth, output_filename, fif, fq);
});

bool PIVOTSBINARY_LCSSORT_REGISTED = registerCallback("pb_lcs_sort", 
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


    typedef ComposibleIndex<PivotsBinary, LCSSort> PB_LCS;

    const auto& fif = [&](){
		auto index = make_unique<PB_LCS>();
        index->initHasher(d, K, M, n, data);
        index->initBucketer(K);
        index->build(n, data);
        return index;
    };

	const auto& fq = [&](PB_LCS &index, int k, const float* queryi, MinK_List* list) {
        int nCandidates = k+sqrt(2*K*M)+K/d*10;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, d, query, ground_truth, output_filename, fif, fq);
});

bool SRP_LCSSORT_REGISTED = registerCallback("srp_lcs_sort", 
		"n qn d K dataset_filename queryset_filename ground_truth_filename output_filename", [&](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int M = sizeof(SigType)*8;
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");


    typedef ComposibleIndex<SRP, LCSSort> SRP_LCS;

    const auto& fif = [&](){
		auto index = make_unique<SRP_LCS>();
        index->initHasher(d, K, M);
        index->initBucketer(K);
        index->build(n, data);
        return index;
    };

	const auto& fq = [&](SRP_LCS &index, int k, const float* queryi, MinK_List* list) {
        int nCandidates = k+K*M;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, d, query, ground_truth, output_filename, fif, fq);
});