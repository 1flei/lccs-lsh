#include "ground_truth.h"

using namespace std;

bool ground_truth_angle_registed = MyCallbackRegister::registerCallback("ground_truth_angle", 
		"n qn d dataset_filename queryset_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	string output_filename = algAs<string>("output_filename");

	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	ground_truth_angle(n, qn, d, data, query, output_filename.c_str());
});

bool ground_truth_cosine_similarity_registed = MyCallbackRegister::registerCallback("ground_truth_cosine_similarity", 
		"n qn d dataset_filename queryset_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	string output_filename = algAs<string>("output_filename");

	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	ground_truth_cosine_similarity(n, qn, d, data, query, output_filename.c_str());
});

bool ground_truth_l2_registed = MyCallbackRegister::registerCallback("ground_truth_l2", 
		"n qn d dataset_filename queryset_filename output_filename", [](){
	using namespace MyCallbackRegister;
	int n = algAs<int>("n");
	int qn = algAs<int>("qn");
	int d = algAs<int>("d");
	string output_filename = algAs<string>("output_filename");

	const float** data = algAs<const float**>("dataset");
	const float** query = algAs<const float**>("queryset");
	ground_truth_l2(n, qn, d, data, query, output_filename.c_str());
});

// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	// timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------

    FILE *fp = fopen(truth_set, "w");
    if (!fp) {
        printf("Could not create %s.\n", truth_set);
        return 1;
    }

    double duration = MyTimer::measure([&](){
        auto list = make_unique<MaxK_List>(MAXK);
        // MaxK_List *list = new MaxK_List(MAXK);
        fprintf(fp, "%d %d\n", qn, MAXK);
        for (int i = 0; i < qn; ++i) {
            list->reset();
            for (int j = 0; j < n; ++j) {	
                float ip = calc_inner_product(d, data[j], query[i]);
                list->insert(ip, j + 1);
            }

            for (int j = 0; j < MAXK; ++j) {
                fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
            }
            fprintf(fp, "\n");
        }
    });
    fclose(fp);

	printf("Ground Truth: %f Seconds\n\n", duration);

	return 0;
}

// -----------------------------------------------------------------------------
int ground_truth_l2(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

    double truth_time = MyTimer::measure([&](){
        auto list = make_unique<MinK_List>(MAXK);
        fprintf(fp, "%d %d\n", qn, MAXK);
        for (int i = 0; i < qn; ++i) {
            list->reset();
            for (int j = 0; j < n; ++j) {
                float ip = calc_l2_dist(d, data[j], query[i]);
                list->insert(ip, j + 1);
            }

            for (int j = 0; j < MAXK; ++j) {
                fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
            }
            fprintf(fp, "\n");
        }
    });
    fclose(fp);

	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}


// -----------------------------------------------------------------------------
int ground_truth_for_weighted_space(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const float **weight,				// query set
	const char  *truth_set)				// address of truth set
{
	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

    double truth_time = MyTimer::measure([&](){
        auto list = make_unique<MinK_List>(MAXK);
        fprintf(fp, "%d %d\n", qn, MAXK);
        for (int i = 0; i < qn; ++i) {
            list->reset();
            for (int j = 0; j < n; ++j) {
                float ip = calc_weighted_dist2(d, weight[i], data[j], query[i]);
                list->insert(ip, j + 1);
            }

            for (int j = 0; j < MAXK; ++j) {
                fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
            }
            fprintf(fp, "\n");
        }
    });
    fclose(fp);

	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}
// -----------------------------------------------------------------------------
int ground_truth_furthest(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

    double truth_time = MyTimer::measure([&](){
        auto list = make_unique<MaxK_List>(MAXK);
        fprintf(fp, "%d %d\n", qn, MAXK);
        for (int i = 0; i < qn; ++i) {
            list->reset();
            for (int j = 0; j < n; ++j) {
                float ip = calc_l2_dist(d, data[j], query[i]);
                list->insert(ip, j + 1);
            }

            for (int j = 0; j < MAXK; ++j) {
                fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
            }
            fprintf(fp, "\n");
        }
    });
	fclose(fp);

	printf("Ground Truth: %f Seconds\n\n", truth_time);
	return 0;
}
// -----------------------------------------------------------------------------
int ground_truth_angle(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

    double truth_time = MyTimer::measure([&](){
        auto list = make_unique<MinK_List>(MAXK);
        fprintf(fp, "%d %d\n", qn, MAXK);
        for (int i = 0; i < qn; ++i) {
            list->reset();
            for (int j = 0; j < n; ++j) {
                float ip = calc_angle(d, data[j], query[i]);
                list->insert(ip, j + 1);
            }

            for (int j = 0; j < MAXK; ++j) {
                fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
            }
            fprintf(fp, "\n");
        }
    });
	fclose(fp);

	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}
// -----------------------------------------------------------------------------
int ground_truth_cosine_similarity(					// find the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

    double truth_time = MyTimer::measure([&](){
        auto list = make_unique<MaxK_List>(MAXK);
        fprintf(fp, "%d %d\n", qn, MAXK);
        for (int i = 0; i < qn; ++i) {
            list->reset();
            for (int j = 0; j < n; ++j) {
                float ip = calc_cosangle(d, data[j], query[i]);
                list->insert(ip, j + 1);
            }

            for (int j = 0; j < MAXK; ++j) {
                fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
            }
            fprintf(fp, "\n");
        }
    });
	fclose(fp);

	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}