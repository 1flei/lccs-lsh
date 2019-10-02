#include "composible_index.h"


using namespace std;
using namespace MyCallbackRegister;

bool SRP_LCSSORT_INT_REGISTED = registerCallback("lcsb",
		"n qn d L r step dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
    using namespace mylccs;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int L = algAs<int>("L");
    double r = algAs<double>("r");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");

    // int checked_candidate = algAs<int>("checked_candidate");
    // double alpha = algAs<double>("alpha");
    int step = algAs<int>("step");

    // using LCSIndex = LCCS_LSH_REORDER;
    using LCSIndex = LCCS_SORT_INT;
    typedef ComposibleIndex<E2Eigen, LCSIndex, int32_t> SRP_LCS;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    double rho = log(L)/log(n);
    // int step = std::max(1, int(exp(3.2-4*rho)));
    // int step = 5;

    printf("rho=%f, step=%d\n", rho, step);

    const auto& fif = [&](){
		auto index = make_unique<SRP_LCS>();
        index->initHasher(d, L*step, r);
        index->initBucketer(L, step);
        index->build(n, data);
        return index;
    };

    fprintf(fp.get(), "lcsb  r=%f, L=%d\n", r, L);
    std::vector<int> checked_candidates = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512};
    // std::vector<int> checked_candidates = {8};
    const auto& fq = [&](SRP_LCS &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        // int nCandidates = k + checked_candidate;
        int nCandidates = checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_l2_dist(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});


bool SRP_LCSSORT_INT_REORDER_REGISTED = registerCallback("lcsb_reorder",
		"n qn d L r step dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
    using namespace mylccs;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int L = algAs<int>("L");
    double r = algAs<double>("r");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");

    // int checked_candidate = algAs<int>("checked_candidate");
    // double alpha = algAs<double>("alpha");
    int step = algAs<int>("step");

    using LCSIndex = LCCS_LSH_REORDER;
    typedef ComposibleIndex<E2Eigen, LCSIndex, int32_t> SRP_LCS;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    double rho = log(L)/log(n);
    // int step = std::max(1, int(exp(3.2-4*rho)));
    // int step = 5;

    printf("rho=%f, step=%d\n", rho, step);

    const auto& fif = [&](){
		auto index = make_unique<SRP_LCS>();
        index->initHasher(d, L*step, r);
        index->initBucketer(L, step);
        index->build(n, data);
        return index;
    };

    fprintf(fp.get(), "lcsb_reorder  r=%f, L=%d\n", r, L);
    std::vector<int> checked_candidates = {1, 2, 4, 8, 16, 32, 64};
    // std::vector<int> checked_candidates = {8};
    const auto& fq = [&](SRP_LCS &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        // int nCandidates = k + checked_candidate;
        int nCandidates = checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_l2_dist(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});

bool E2LSH_REGISTED = registerCallback("e2lsh",
		"n qn d K L r dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int L = algAs<int>("L");
    double r = algAs<double>("r");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef ComposibleIndex<E2Eigen, KLBucketing> Index;
    typedef ComposibleIndex<E2Eigen, KLBucketingSimpleHasher> Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>();
        index->initHasher(d, K*L, r);
        index->initBucketer(n, L, K);
        index->inc_build(n, data);
        return index;
    };

    fprintf(fp.get(), "e2lsh  r=%f, K=%d L=%d\n", r, K, L);
    std::vector<int> checked_candidates = {16, 64, 256, 1024, 4096, 16384, 65536, 262144};
    // std::vector<int> checked_candidates;
    // for(int check_k = 64; check_k < n/10; check_k*=4){
    //     checked_candidates.push_back(check_k);
    // }
    const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        int nCandidates = k + checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_l2_dist(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});

// bool MPLSH_REGISTED = registerCallback("mplsh",
// 		"n qn d K L r dataset_filename queryset_filename ground_truth_filename output_filename", [](){
// 	using namespace MyCallbackRegister;
// 	int n = algAs<int>("n");
// 	int qn = algAs<int>("qn");
// 	int d = algAs<int>("d");
// 	int K = algAs<int>("K");
//     int L = algAs<int>("L");
//     double r = algAs<double>("r");
// 	const float** data = algAs<const float**>("dataset");
// 	const float** query = algAs<const float**>("queryset");
// 	const Result** ground_truth = algAs<const Result**>("ground_truth");
// 	string output_filename = algAs<string>("output_filename");
//     // int checked_candidate = algAs<int>("checked_candidate");


//     typedef ComposibleIndex<E2Eigen, MYMPLSH> Index;
//     std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

//     const auto& fif = [&](){
// 		auto index = make_unique<Index>();
//         index->initHasher(d, K*L, r);
//         index->initBucketer(n, L, K);
//         index->inc_build(n, data);
//         return index;
//     };

//     fprintf(fp.get(), "mplsh  r=%f, K=%d L=%d\n", r, K, L);
//     // std::vector<int> checked_candidates = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
//     std::vector<int> checked_candidates;
//     for(int check_k = 1; check_k <= 64; check_k*=2){
//         checked_candidates.push_back(check_k);
//     }
//     const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
//         // int nCandidates = k+K*M;
//         int nCandidates = k + checked_candidate;
//         const auto& f = [&](int idx){
//             float angle = calc_l2_dist(d, data[idx], queryi);
//             list->insert(angle, idx + 1);
//         };
//         index.query(nCandidates, queryi, f);
//     };

//     benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
// });

bool C2LSH_REGISTED = registerCallback("c2lsh",
		"n qn d L r dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	// int K = algAs<int>("K");
    int L = algAs<int>("L");
    double r = algAs<double>("r");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    typedef ComposibleIndex<E2Eigen, C2Bucketing> Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>();
        index->initHasher(d, L, r);
        index->initBucketer(n, L);
        index->inc_build(n, data);
        
        return index;
    };

    fprintf(fp.get(), "c2lsh  r=%f, L=%d\n", r, L);
    std::vector<int> checked_candidates = {64, 128, 256, 512, 1024, 2048};

    const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        int nCandidates = k + checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_l2_dist(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});


bool SRP_E2_REGISTERED = registerCallback("srp_e2",
		"n qn d K L dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int L = algAs<int>("L");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef ComposibleIndex<E2Eigen, KLBucketing> Index;
    typedef ComposibleIndex<SRPCompact, KLBucketingCompact, uint64_t> Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>();
        index->initHasher(d, K*L);
        index->initBucketer(n, L, K);
        index->inc_build(n, data);
        return index;
    };

    fprintf(fp.get(), "srp_e2, K=%d L=%d\n", K, L);
    // std::vector<int> checked_candidates = {16, 64, 256, 1024, 4096, 16384, 65536, 262144};
    std::vector<int> checked_candidates;
    for(int check_k = 64; check_k < n/10; check_k*=2){
        checked_candidates.push_back(check_k);
    }
    const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        int nCandidates = k + checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});


bool POLYTOPE_E2_REGISTERED = registerCallback("polytope_e2",
		"n qn d K L dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	int K = algAs<int>("K");
    int L = algAs<int>("L");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef ComposibleIndex<E2Eigen, KLBucketing> Index;
    typedef ComposibleIndex<PolytopeHasher, KLBucketing> Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>();
        index->initHasher(d, K*L);
        index->initBucketer(n, L, K);
        index->inc_build(n, data);
        return index;
    };

    fprintf(fp.get(), "polytope_e2, K=%d L=%d\n", K, L);
    // std::vector<int> checked_candidates = {16, 64, 256, 1024, 4096, 16384, 65536, 262144};
    std::vector<int> checked_candidates;
    for(int check_k = 1; check_k < n/10; check_k*=2){
        checked_candidates.push_back(check_k);
    }
    const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        int nCandidates = k + checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});

bool POLYTOPE_C2_REGISTERED = registerCallback("polytope_c2",
		"n qn d L dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
    int L = algAs<int>("L");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef ComposibleIndex<E2Eigen, KLBucketing> Index;
    typedef ComposibleIndex<PolytopeHasher, C2Bucketing> Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>();
        index->initHasher(d, L);
        index->initBucketer(n, L);
        index->inc_build(n, data);
        return index;
    };

    fprintf(fp.get(), "polytope_e2, L=%d\n", L);
    std::vector<int> checked_candidates;
    for(int check_k = 1; check_k < n/10; check_k*=2){
        checked_candidates.push_back(check_k);
    }
    const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        int nCandidates = k + checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});

bool POLYTOPE_LCCS_REGISTERED = registerCallback("polytope_lccs",
		"n qn d L step dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
    int L = algAs<int>("L");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");
    int step = algAs<int>("step");

    // for(int i=1;i<=qn;i*=10){
    //     printf("q%d=\n", i);
    //     printVec(query[i-1], 300);
    // }

    // typedef ComposibleIndex<E2Eigen, KLBucketing> Index;
    typedef ComposibleIndex<PolytopeHasher, mylccs::LCCS_SORT_INT> Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>();
        index->initHasher(d, L);
        index->initBucketer(L, step);
        index->build(n, data);
        return index;
    };


    fprintf(fp.get(), "polytope_lccs, L=%d\n", L);
    std::vector<int> checked_candidates;
    for(int check_k = 1; check_k <= 1024; check_k*=2){
        checked_candidates.push_back(check_k);
    }
    // std::vector<int> checked_candidates = {1};
    const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        int nCandidates = checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});

bool SRP_SCAN_REGISTERED = registerCallback("srp_scan",
		"n qn d L dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
    int L = algAs<int>("L");
	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	const Result** ground_truth = algAs<const Result**>("ground_truth");
	string output_filename = algAs<string>("output_filename");
    // int checked_candidate = algAs<int>("checked_candidate");


    // typedef ComposibleIndex<E2Eigen, KLBucketing> Index;
    typedef ComposibleIndex<SRPCompact, HammingLinearScan, uint64_t> Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    const auto& fif = [&](){
		auto index = make_unique<Index>();
        index->initHasher(d, L);
        index->initBucketer(L);
        index->build(n, data);
        return index;
    };

    fprintf(fp.get(), "srp_scan, L=%d\n", L);
    std::vector<int> checked_candidates = {16};
    // std::vector<int> checked_candidates;
    // for(int check_k = 64; check_k < n/10; check_k*=4){
    //     checked_candidates.push_back(check_k);
    // }
    const auto& fq = [&](Index &index, int k, int checked_candidate, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        int nCandidates = k + checked_candidate;
        const auto& f = [&](int idx){
            float angle = calc_angle(d, data[idx], queryi);
            printf("     %d:%f\n", idx, angle);
            list->insert(angle, idx + 1);
        };
        index.query(nCandidates, queryi, f);
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});