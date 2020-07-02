#include "qalsh.h"

using namespace std;
using namespace MyCallbackRegister;


// -----------------------------------------------------------------------------
// bad implementation. Required by QALSH but including util.cc will conflict with my util.cc for some methods


// -----------------------------------------------------------------------------
float calc_lp_dist(					// calc L_{p} norm
	int   dim,							// dimension
	float p,							// the p value of Lp norm, p in (0,2]
	float threshold,					// threshold
	const float *vec1,					// 1st point
	const float *vec2)					// 2nd point
{
	if (fabs(p - 2.0f) < FLOATZERO) {
		return sqrt(calc_l2_sqr(dim, vec1, vec2));
	}
	else if (fabs(p - 1.0f) < FLOATZERO) {
		return calc_l1_dist(dim, vec1, vec2);
	}
	else if (fabs(p - 0.f) < FLOATZERO) {
		return calc_l0_dist(dim, vec1, vec2);
	}
	else {
		return calc_lp_dist_p(dim, p, vec1, vec2);
	}
}


bool QALSH_REGISTERED = registerCallback("qalsh",
		"n qn d L M c dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	int n = argAs<int>("n");
	int qn = argAs<int>("qn");
	int d = argAs<int>("d");

	const float** data = argAs<const float**>("dataset");
	const float** query = argAs<const float**>("queryset");
	const Result** ground_truth = argAs<const Result**>("ground_truth");
	std::string output_filename = argAs<std::string>("output_filename");
    // int checked_candidate = argAs<int>("checked_candidate");

    // typedef MPLSH Index;
    typedef QALSH_PLUS Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);

    //params confirmed by authors
    int leaf = 20000;
    int L = argAs<int>("L");
    int M = argAs<int>("M");
    double p = 2.;
    double zeta = 0.;
    double ratio = argAs<double>("c");

    const auto& fif = [&](){
        auto index = make_unique<Index>(n, d, leaf, L, M, p, zeta, ratio, data);
        return index;
    };

    fprintf(fp.get(), "qalsh1000 L=%d M=%d\n", L, M);
    // std::vector<int> ts = {1, 4, 16, 64, 256, 2048, 8192};
    // std::vector<int> ts = {64};
    // std::vector<int> nbs = {2, 4, 8, 16, 32, 64};
    std::vector<int> nbs = {2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 48, 64};

    const auto& fq = [&](Index &index, int k, int nb, const float* queryi, MinK_List* list) {
        index.knn(k, nb, queryi, list);
        // printf("----------------\n");
        // for(int i=0;i<list->getk();i++){
        //   printf("(%f, %d)\n", list->ith_key(i), list->ith_id(i));
        // }
    };

    benchmarkMinklist(qn, query, ground_truth, nbs, fp.get(), fif, fq);
});