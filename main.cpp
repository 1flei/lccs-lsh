#include <boost/program_options.hpp>
#include <iostream>
#include <string> 
#include "register.h"
#include "myndarray.h"
#include "util.h"
#include "myTimer.h"

using namespace std;
using namespace boost::program_options;

//some initilization and registration will be done before main function
//can check register.h/cpp and composible_index.cpp for example
int main(int argc, char **argv)
{
    int n, qn, d;
	string datasetFilename, queryFilename, weightFilename, groundtruthFilename, outputFilename;

	// srand(time(NULL));
	srand(GLOBAL_SEED);

    // Declare the supported options.
    options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("algorithm_name,A", value<std::string>(), "the name of algorithm to run")

		("n,n", value(&n), "the number of data points")
		("d,d", value(&d), "the dimension of data")
		("qn,q", value(&qn), "the number of query points")
		
        ("dataset_filename,D", value(&datasetFilename), "path to dataset filename")
		("queryset_filename,Q", value(&queryFilename), "path to query filename")
		("ground_truth_filename,G", value(&groundtruthFilename), "path to ground truth filename")
		("output_filename,O", value(&outputFilename), "output folder path (with / at the end) or output filename")

        ("K,K", value<int>(), "parameter used for some algorithms")
        ("M,M", value<int>(), "parameter used for some algorithms")
        ("L,L", value<int>(), "parameter used for some algorithms")
        ("m,m", value<int>(), "parameter used for some algorithms")
        ("p,p", value<double>(), "extra probes multiplication of mp_lccs")
        ("step", value<int>(), "parameter used for some algorithms")
        ("r,r", value<double>(), "parameter used for some algorithms")
        ("c,c", value<double>(), "parameter used for some algorithms")

		("nHashBits", value<int>(), "parameter used for FALCONN")
        // ("cp_dim", value<int>(), "parameter for Crosspolytope Hasher")
		("bit_packed", "whether to use bit_packed hash table for FALCONN")
        ("cnt_threshold", value<int>(), "parameter used for C2LSH")

		("binary_input", "read from binary input")
		("normalized", "whether to normalize the dataset")

		("alpha", value<double>(), "parameter used for some algorithms")

		("srs_targeted_dimension", value<int>(), "parameter used for SRS")
        ("p_thres", value<double>(), "parameter used for SRS")
		
		("checked_candidate", value<int>()->default_value(100), "the number of candidates to verify for each algorithm")
    ;

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);  

    if (vm.count("help")) {
        cout << desc << "\n";

        cout << "available algorithms (and required parameters) are:" << endl;
        MyCallbackRegister::showAllRegisteredCallback();
        return 1;
    }

    
	// -------------------------------------------------------------------------
	//  read whatever needed
	// -------------------------------------------------------------------------
	typedef NDArray<2, Scalar> F2DArray;

	unique_ptr<F2DArray > dataArr, queryArr;
	unique_ptr<NDArray<2, Result> > resultArr;

	float **data=nullptr;
	float **query=nullptr;
	Result **results=nullptr;

	if(vm.count("binary_input")){
		if(datasetFilename!=""){
			dataArr = unique_ptr<F2DArray>(new F2DArray({(size_t)n, (size_t)d}));
			data = dataArr->to_ptr();
			if (read_data_binary(n, d, datasetFilename.c_str(), data) == 1) {
				printf("Reading dataset error!\n");
				return 1;
			}
		}

		if(queryFilename!=""){
			queryArr = unique_ptr<F2DArray>(new F2DArray({(size_t)qn, (size_t)d}));
			query = queryArr->to_ptr();
			if (read_data_binary(qn, d, queryFilename.c_str(), query) == 1) {
				printf("Reading query set error!\n");
				return 1;
			}
		}
	} else{
		if(datasetFilename!=""){
			dataArr = unique_ptr<F2DArray>(new F2DArray({(size_t)n, (size_t)d}));
			data = dataArr->to_ptr();
			if (read_data(n, d, datasetFilename.c_str(), data) == 1) {
				printf("Reading dataset error!\n");
				return 1;
			}
		}

		if(queryFilename!=""){
			queryArr = unique_ptr<F2DArray>(new F2DArray({(size_t)qn, (size_t)d}));
			query = queryArr->to_ptr();
			if (read_data(qn, d, queryFilename.c_str(), query) == 1) {
				printf("Reading query set error!\n");
				return 1;
			}
		}
	}

	if(groundtruthFilename!=""){
		resultArr = unique_ptr<NDArray<2, Result>>(new NDArray<2, Result>({(size_t)qn, (size_t)MAXK}));
		results = resultArr->to_ptr();
		if (read_ground_truth(qn, groundtruthFilename.c_str(), results) == 1) {
			printf("Reading Truth Set Error!\n");
			return 1;
		}
	}
    cout << "finishing reading data, query and ground truth" << endl;

	if(vm.count("normalized")){
		for(int i=0;i<n; i++){
			normalize(d, data[i]);
		}
		for(int i=0;i<qn; i++){
			normalize(d, query[i]);
		}
	}

    using namespace MyCallbackRegister;
    if (!vm.count("algorithm_name")) {
        cout << "algorithm_name required" << endl;

        cout << "see --help for more information" << endl;
    } else{
        // cout << vm["algorithm_name"].as<string>() << endl;
        setvm(vm);
        addArg("dataset", (const float**)data);
        addArg("queryset", (const float**)query);
        addArg("ground_truth", (const Result**)results);
        
        run(vm["algorithm_name"].as<string>());
    }

	MyTimer::printAll();
    return 0;
}