#include "headers.h"

// -----------------------------------------------------------------------------
void usage() 						// display the usage of this package
{
	printf("\n"
		"-------------------------------------------------------------------\n"
		" Usage of the package for Maximum Cosine Similarity Search (MCSS)  \n"
		"-------------------------------------------------------------------\n"
		"    -alg  {integer}  options of algorithms (0 - 13)\n"
		"    -n    {integer}  cardinality of the dataset\n"
		"    -qn   {integer}  number of queries\n"
		"    -d    {integer}  dimensionality of the dataset\n"
		"    -m    {integer}  #hash func for m-SRP and LCSB\n"
		"    -k    {integer}  #concatenations for (k,l)-SRP\n"
		"    -l    {integer}  #hash tables for (k,l)-SRP\n"
		"    -ds   {string}   address of the data  set\n"
		"    -qs   {string}   address of the query set\n"
		"    -ts   {string}   address of the truth set\n"
		"    -of   {string}   output folder\n"
		"\n"
		"-------------------------------------------------------------------\n"
		" The options of algorithms are:\n"
		"-------------------------------------------------------------------\n"
		"    0  - Ground-Truth\n"
		"         Parameters: -alg 0 -n -qn -d -ds -qs -ts\n"
		"\n"
		"    1  - MCSS by Linear Scan\n"
		"         Parameters: -alg 1 -n -qn -d -ds -qs -ts -of\n"
		"\n"
		"    2  - MCSS by m-SRP\n"
		"         Parameters: -alg 2 -n -qn -d -m -ds -qs -ts -of\n"
		"\n"
		"    3  - MCSS by (k,l)-SRP\n"
		"         Parameters: -alg 2 -n -qn -d -k -l -ds -qs -ts -of\n"
		"\n"
		"    4  - MCSS by Longest Common Substring Bucketing (LCSB)\n"
		"         Parameters: -alg 4 -n -qn -d -m -ds -qs -ts -of\n"
		"\n"
		"-------------------------------------------------------------------\n"
		" Authors: Qiang Huang (huangq2011@gmail.com)                       \n"
		"          Yifan Lei   (yfleiii@gmail.com)                          \n"
		"-------------------------------------------------------------------\n"
		"\n\n\n");
}

// -----------------------------------------------------------------------------
int main(int nargs, char **args)
{
	srand(666);						// set the random seed
	//usage();

	int   alg = -1;					// which algorithm?
	int   n   = -1;					// cardinality
	int   qn  = -1;					// query number
	int   d   = -1;					// dimensionality
	int   m   = -1;					// #hash func for m-SRP and LCSB
	int   k   = -1;					// #concatenations for (k,l)-SRP
	int   l   = -1;					// #hash tables for (k,l)-SRP
	
	char  data_set[200];			// address of data set
	char  query_set[200];			// address of query set
	char  truth_set[200];			// address of ground truth file
	char  output_folder[200];		// output folder

	bool  failed = false;
	int   cnt = 1;
	
	while (cnt < nargs && !failed) {
		if (strcmp(args[cnt], "-alg") == 0) {
			alg = atoi(args[++cnt]);
			printf("alg           = %d\n", alg);
			if (alg < 0 || alg > 11) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-n") == 0) {
			n = atoi(args[++cnt]);
			printf("n             = %d\n", n);
			if (n <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-d") == 0) {
			d = atoi(args[++cnt]);
			printf("d             = %d\n", d);
			if (d <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-qn") == 0) {
			qn = atoi(args[++cnt]);
			printf("qn            = %d\n", qn);
			if (qn <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-m") == 0) {
			m = atoi(args[++cnt]);
			printf("m             = %d\n", m);
			if (m <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-k") == 0) {
			k = atoi(args[++cnt]);
			printf("k             = %d\n", k);
			if (k <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-l") == 0) {
			l = atoi(args[++cnt]);
			printf("l             = %d\n", l);
			if (l <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-ds") == 0) {
			strncpy(data_set, args[++cnt], sizeof(data_set));
			printf("data_set      = %s\n", data_set);
		}
		else if (strcmp(args[cnt], "-qs") == 0) {
			strncpy(query_set, args[++cnt], sizeof(query_set));
			printf("query_set     = %s\n", query_set);
		}
		else if (strcmp(args[cnt], "-ts") == 0) {
			strncpy(truth_set, args[++cnt], sizeof(truth_set));
			printf("truth_set     = %s\n", truth_set);
		}
		else if (strcmp(args[cnt], "-of") == 0) {
			strncpy(output_folder, args[++cnt], sizeof(output_folder));
			printf("output_folder = %s\n", output_folder);

			int len = (int) strlen(output_folder);
			if (output_folder[len - 1] != '/') {
				output_folder[len] = '/';
				output_folder[len + 1] = '\0';
			}
			create_dir(output_folder);
		}
		else {
			failed = true;
			usage();
			break;
		}
		cnt++;
	}
	printf("\n");

	// -------------------------------------------------------------------------
	//  read data set and query set
	// -------------------------------------------------------------------------
	timeval start_time, end_time;

	gettimeofday(&start_time, NULL);
	float** data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading dataset error!\n");
		return 1;
	}

	float** query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading query set error!\n");
		return 1;
	}
	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data and Query: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  methods
	// -------------------------------------------------------------------------
	switch (alg) {
	case 0:
		ground_truth(n, qn, d, (const float **) data, (const float **) query, 
			truth_set);
		break;
	case 1:
		linear_scan(n, qn, d, (const float **) data, (const float **) query,
			truth_set, output_folder);
		break;
	case 2:
		m_srp(n, qn, d, m, (const float **) data, (const float **) query, 
			truth_set, output_folder);
		break;
	case 3:
		kl_srp(n, qn, d, k, l, (const float **) data, (const float **) query, 
			truth_set, output_folder);
		break;
	case 4:
		lcsb(n, qn, d, m, (const float **) data, (const float **) query, 
			truth_set, output_folder);
		break;
	default:
		printf("Parameters error!\n");
		usage();
		break;
	}
	
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data  = NULL;

	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;

	// uint64_t a = (uint64_t) 1 << 60;
	// printf("a = %llu\n", a);

	// uint64_t b = (a >> 16) << 16;
	// printf("b = %llu\n", b);

	// int c = 0 % 64;
	// printf("c = %d\n", c);

	// unordered_set<int> s;

	// if (s.find(3) == s.end()) {
	// 	printf("not found\n");
	// 	s.insert(3);
	// }

	// // s.insert(4);
	// for (auto it = s.begin(); it != s.end(); ++it) {
	// 	printf("%d ", *it);
	// }
	// printf("\n");

	// if (s.find(3) != s.end()) {
	// 	printf("we find 3\n");
	// }


	// map<vector<int>, vector<int> > ht;

	// vector<int> a;
	// a.push_back(1); a.push_back(2); a.push_back(3);

	// ht.insert()
	// ht[a].push_back(905);
	// ht[a].push_back(901);
	// ht[a].push_back(908);

	// vector<int> b;
	// b.push_back(1); b.push_back(3);
	// ht[a].push_back(805);
	// ht[a].push_back(801);
	// ht[a].push_back(808);

	// for(auto it = ht.begin(); it != ht.end(); it++){
    //     vector<int> key = it->first;
	// 	vector<int> val = it->second;

	// 	printf("key:");
	// 	for (int i = 0; i < key.size(); ++i) printf(" %d", key[i]);
	// 	printf(", val:");
	// 	for (int i = 0; i < val.size(); ++i) printf(" %d", val[i]);
    // }

	// Pri_List *list = new Pri_List(5);

	// list->insert(98, 5);
	// list->insert(96, 1);
	// list->insert(88, 2);
	// list->insert(82, 3);
	// list->insert(60, 4);
	// for (int i = 0; i < 5; ++i) {
	// 	printf("(%d, %d) -> ", list->ith_key(i), list->ith_id(i));
	// }
	// printf("\n");

	// list->insert(92, 3);
	// for (int i = 0; i < 5; ++i) {
	// 	printf("(%d, %d) -> ", list->ith_key(i), list->ith_id(i));
	// }
	// printf("\n");

	// list->insert(85, 2);
	// for (int i = 0; i < 5; ++i) {
	// 	printf("(%d, %d) -> ", list->ith_key(i), list->ith_id(i));
	// }
	// printf("\n");

	// list->insert(98, 6);
	// for (int i = 0; i < 5; ++i) {
	// 	printf("(%d, %d) -> ", list->ith_key(i), list->ith_id(i));
	// }
	// printf("\n");

	// list->insert(99, 7);
	// for (int i = 0; i < 5; ++i) {
	// 	printf("(%d, %d) -> ", list->ith_key(i), list->ith_id(i));
	// }
	// printf("\n");

	// list->insert(40, 8);
	// for (int i = 0; i < 5; ++i) {
	// 	printf("(%d, %d) -> ", list->ith_key(i), list->ith_id(i));
	// }
	// printf("\n");

	return 0;
}


