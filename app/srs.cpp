#include "srs.h"

using namespace std;
using namespace MyCallbackRegister;


double cal_thres(double c, double p_thres, int m) {
  if (p_thres >= 1) {
    return -1;
  }
  boost::math::chi_squared chi(m);
  return boost::math::quantile(chi, p_thres) / c / c;
}


// bool SRS_REGISTERED = registerCallback("srs",
// 		"n qn d c r srs_targeted_dimension dataset_filename queryset_filename ground_truth_filename output_filename", [](){
// 	int n = argAs<int>("n");
// 	int qn = argAs<int>("qn");
// 	int d = argAs<int>("d");
//     int srs_targeted_dimension = argAs<int>("srs_targeted_dimension");
//     double c = argAs<double>("c");
//     double p_thres = argAs<double>("p_thres");

// 	const float** data = argAs<const float**>("dataset");
// 	const float** query = argAs<const float**>("queryset");
// 	const Result** ground_truth = argAs<const Result**>("ground_truth");
// 	std::string output_filename = argAs<std::string>("output_filename");
//     // int checked_candidate = argAs<int>("checked_candidate");

//     // typedef MPLSH Index;
//     typedef SRS_In_Memory<float> Index;
//     std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);


//     const auto& fif = [&](){

// 		auto index = make_unique<Index>("./");
//         // index->build_index_from_data(n, d, srs_targeted_dimension, (float**)data);
//         // index->restore_index();
//         index->build_index_in_memory(n, d, srs_targeted_dimension, (float**)data);
//         return index;
//     };

//     fprintf(fp.get(), "srs  c=%f, p_thres=%f, targeted_dimension=%d\n", c, p_thres, srs_targeted_dimension);
//     // std::vector<int> ts = {1, 4, 16, 64, 256, 2048, 8192};
//     std::vector<int> ts = {64};

//     const auto& fq = [&](Index &index, int k, int t, const float* queryi, MinK_List* list) {
//         // int nCandidates = k+K*M;
//         typedef typename Accumulator<float>::Type ResultType;
//         std::vector<res_pair_raw<ResultType> > res;
//         // index.knn_search((float*)queryi, k, t + k - 1, cal_thres(c, p_thres, srs_targeted_dimension), res);
//         index.knn_search((float*)queryi, k, t + k - 1, -1, res);
//         std::cout << "--------------------------" << std::endl;
//         for(auto& resi:res){
//             list->insert(sqrt(resi.dist), resi.id+1);
//             std::cout << sqrt(resi.dist) << ", " << resi.id+1 << std::endl;
//         }
//     };

//     benchmarkMinklist(qn, query, ground_truth, ts, fp.get(), fif, fq);
// });





bool SRS_REGISTERED = registerCallback("srs",
		"n qn d L dataset_filename queryset_filename ground_truth_filename output_filename", [](){
	int n = argAs<int>("n");
	int qn = argAs<int>("qn");
	int d = argAs<int>("d");
    int srs_targeted_dimension = argAs<int>("L");
    // double c = argAs<double>("c");
    // double p_thres = argAs<double>("p_thres");

	const float** data = argAs<const float**>("dataset");
	const float** query = argAs<const float**>("queryset");
	const Result** ground_truth = argAs<const Result**>("ground_truth");
	std::string output_filename = argAs<std::string>("output_filename");
    // int checked_candidate = argAs<int>("checked_candidate");

    // typedef MPLSH Index;
    typedef MYSRS Index;
    std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(output_filename.c_str(), "a+"), &fclose);


    const auto& fif = [&](){
		  auto index = make_unique<Index>(n, d, srs_targeted_dimension);
      // index->build_index_from_data(n, d, srs_targeted_dimension, (float**)data);
      // index->restore_index();
      index->build(data);
      return index;
    };

    fprintf(fp.get(), "srs L=%d\n", srs_targeted_dimension);
    // std::vector<int> ts = {1, 4, 16, 64, 256, 2048, 8192};
    // std::vector<int> ts = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
    std::vector<int> checked_candidates;
    for(int check_k = 64; check_k < n/2; check_k*=2){
        checked_candidates.push_back(check_k);
    }
    // std::vector<int> ts = {64};

    const auto& fq = [&](Index &index, int k, int t, const float* queryi, MinK_List* list) {
        // int nCandidates = k+K*M;
        index.query(t + k - 1, queryi, list);
        // printf("----------------\n");
        // for(int i=0;i<list->getk();i++){
          // printf("(%f, %d)\n", list->ith_key(i), list->ith_id(i));
        // }
    };

    benchmarkMinklist(qn, query, ground_truth, checked_candidates, fp.get(), fif, fq);
});