#include "headers.h"

// -----------------------------------------------------------------------------
int ground_truth(					// find the ground truth results
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
	gettimeofday(&start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	MaxK_List *list = new MaxK_List(MAXK);
	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		for (int j = 0; j < n; ++j) {
			float cosine = calc_cosine(d, data[j], query[i]);
			list->insert(cosine, j + 1);
		}

		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&end_time, NULL);
	float truth_time = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
		start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete list; list = NULL;
	
	return 0;
}


// -----------------------------------------------------------------------------
int linear_scan(					// find top-k mip using linear_scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  mcss by linear scan
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%slinear.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kMCSS[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	printf("Top-k MCSS of Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < max_round; num++) {
		gettimeofday(&start_time, NULL);
		top_k = kMCSS[num];
		MaxK_List* list = new MaxK_List(top_k);

		overall_ratio = 0.0f;
		recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			for (int j = 0; j < n; ++j) {
				float cosine = calc_cosine(d, data[j], query[i]);
				list->insert(cosine, j + 1);
			}
			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			int   size  = list->size();
			for (int j = 0; j < size; ++j) {
				ratio += R[i][j].key_ / list->ith_key(j);
			}
			overall_ratio += ratio / size;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall        = recall / qn;
		runtime       = (runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
			runtime, recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] R; R = NULL;

	return 0;
}


// -----------------------------------------------------------------------------
int m_srp(							// m-srp for MCSS
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// number of projections (concatenation)
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	M_SRP *srp = new M_SRP(n, d, m, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  mcss by m-SRP
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%sm_srp.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	fprintf(fp, "m = %d\n", m);
	fprintf(fp, "indexing time = %.6f\n", indexing_time);

	int kMCSS[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	int kCHECKS[] = {50, 100, 200, 500, 1000, 2000, 5000, 10000};
	int max_check = 8;
	int check_k = -1;
	int check_threshold = (int) floor(0.01f * n);

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	for (int check = 0; check < max_check; ++check) {
		check_k = kCHECKS[check];
		if (check_k > check_threshold) break;
		fprintf(fp, "check_k = %d\n", check_k);

		printf("Top-k MCSS of M_SRP (check_k = %d): \n", check_k);
		printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
		for (int round = 0; round < max_round; round++) {
			gettimeofday(&start_time, NULL);
			top_k = kMCSS[round];
			MaxK_List *list = new MaxK_List(top_k);

			overall_ratio = 0.0f;
			recall = 0.0f;
			for (int i = 0; i < qn; ++i) {
				list->reset();
				srp->mcss(top_k, check_k, query[i], list);
				recall += calc_recall(top_k, (const Result *) R[i], list);
				
				float ratio = 0.0f;
				int   size  = list->size();
				for (int j = 0; j < size; ++j) {
					ratio += R[i][j].key_ / list->ith_key(j);
				}
				overall_ratio += ratio / size;
			}
			delete list; list = NULL;
			gettimeofday(&end_time, NULL);
			runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
				start_time.tv_usec) / 1000000.0f;

			overall_ratio = overall_ratio / qn;
			recall        = recall / qn;
			runtime       = (runtime * 1000.0f) / qn;

			printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
				runtime, recall);
			fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete srp; srp = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int kl_srp(							// (k,l)-srp for MCSS
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   k,							// number of projections (concatenation)
	int   l,							// number of hash tables (repetion)
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	KL_SRP *srp = new KL_SRP(n, d, k, l, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  mcss by (k,l)-SRP
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%skl_srp.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	fprintf(fp, "k = %d, l = %d\n", k, l);
	fprintf(fp, "indexing time = %.6f\n\n", indexing_time);

	int kMCSS[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	int kCHECKS[] = {50, 100, 200, 500, 1000, 2000, 5000, 10000};
	int max_check = 8;
	int check_k = -1;
	int check_threshold = (int) floor(0.01f * n);

	int   count = 0;
	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;
	float hit_rate = -1.0f;

	for (int check = 0; check < max_check; ++check) {
		check_k = kCHECKS[check];
		if (check_k > check_threshold) break;
		fprintf(fp, "check_k = %d\n", check_k);

		printf("Top-k MCSS of KL_SRP (check_k = %d): \n", check_k);
		printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\t\tHit_Rate\n");
		for (int round = 0; round < max_round; round++) {
			gettimeofday(&start_time, NULL);
			top_k = kMCSS[round];
			MaxK_List *list = new MaxK_List(top_k);

			count = 0;
			overall_ratio = 0.0f;
			recall = 0.0f;
			for (int i = 0; i < qn; ++i) {
				list->reset();
				srp->mcss(top_k, check_k, query[i], list);
				recall += calc_recall(top_k, (const Result *) R[i], list);
				
				float ratio = 0.0f;
				int   size  = list->size();
				if (size > 0) {
					for (int j = 0; j < size; ++j) {
						ratio += R[i][j].key_ / list->ith_key(j);
					}
					overall_ratio += ratio / size;
					count++;
				}
			}
			delete list; list = NULL;
			gettimeofday(&end_time, NULL);
			runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
				start_time.tv_usec) / 1000000.0f;

			overall_ratio = overall_ratio / count;
			recall        = recall / count;
			runtime       = (runtime * 1000.0f) / qn;
			hit_rate      = (count * 100.0f) / qn;

			printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\t\t%.2f%%\n", top_k, 
				overall_ratio, runtime, recall, hit_rate);
			fprintf(fp, "%d\t%f\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, 
				recall, hit_rate);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete srp; srp = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int lcsb(							// lcsb for MCSS
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   m,							// number of projections (concatenation)
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set,				// address of truth set
	const char  *output_folder) 		// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read the ground truth file
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	LCSB *lcsb = new LCSB(n, d, m, data);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);
	
	// -------------------------------------------------------------------------
	//  mcss by LCSB
	// -------------------------------------------------------------------------	
	char output_set[200];
	sprintf(output_set, "%slcsb.out", output_folder);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}
	fprintf(fp, "m = %d\n", m);
	fprintf(fp, "indexing time = %.6f\n\n", indexing_time);

	int kMCSS[] = { 1, 2, 5, 10 };
	int max_round = 4;
	int top_k = -1;

	int kCHECKS[] = {50, 100, 200, 500, 1000, 2000, 5000, 10000};
	int max_check = 8;
	int check_k = -1;
	int check_threshold = (int) floor(0.01f * n);

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;

	for (int check = 0; check < max_check; ++check) {
		check_k = kCHECKS[check];
		if (check_k > check_threshold) break;
		fprintf(fp, "check_k = %d\n", check_k);

		printf("Top-k MCSS of LCSB (check_k = %d): \n", check_k);
		printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
		for (int round = 0; round < max_round; round++) {
			gettimeofday(&start_time, NULL);
			top_k = kMCSS[round];
			MaxK_List *list = new MaxK_List(top_k);

			overall_ratio = 0.0f;
			recall = 0.0f;
			for (int i = 0; i < qn; ++i) {
				list->reset();
				lcsb->mcss(top_k, check_k, query[i], list);
				recall += calc_recall(top_k, (const Result *) R[i], list);
				
				float ratio = 0.0f;
				int   size  = list->size();
				for (int j = 0; j < size; ++j) {
					ratio += R[i][j].key_ / list->ith_key(j);
				}
				overall_ratio += ratio / size;
			}
			delete list; list = NULL;
			gettimeofday(&end_time, NULL);
			runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
				start_time.tv_usec) / 1000000.0f;

			overall_ratio = overall_ratio / qn;
			recall        = recall / qn;
			runtime       = (runtime * 1000.0f) / qn;

			printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio, 
				runtime, recall);
			fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, overall_ratio, runtime, recall);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lcsb; lcsb = NULL;
	delete[] R; R = NULL;

	return 0;
}

