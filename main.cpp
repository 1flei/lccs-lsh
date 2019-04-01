#include <boost/program_options.hpp>
#include <iostream>
#include <string> 
#include "register.h"
#include "myndarray.h"
#include "util.h"

using namespace std;
using namespace boost::program_options;


int main(int argc, char **argv)
{
    int n, qn, d;
	string datasetFilename, queryFilename, weightFilename, groundtruthFilename, outputFilename;

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
        
        ("msg", value<string>(), "msg")
    ;

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);  

    if (vm.count("help")) {
        cout << desc << "\n";

        cout << "available algorithms are:" << endl;
        MyCallbackRegister::showAllRegisteredCallback();
        return 1;
    }

    
	// -------------------------------------------------------------------------
	//  read whatever needed
	// -------------------------------------------------------------------------
	typedef NDArray<2, Scalar> F2DArray;

	unique_ptr<F2DArray > dataArr, queryArr, weightArr;
	unique_ptr<NDArray<2, Result> > resultArr;

	float **data=nullptr;
	float **query=nullptr;
	float **weight=nullptr;
	Result **results=nullptr;

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

	if(weightFilename!=""){
		weightArr = unique_ptr<F2DArray>(new F2DArray({(size_t)qn, (size_t)d}));
		weight = weightArr->to_ptr();
		if (read_data(qn, d, weightFilename.c_str(), weight) == 1) {
			printf("Reading query weight error!\n");
			printf("Will not use weight then!\n");
			return 1;
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



    using namespace MyCallbackRegister;
    if (!vm.count("algorithm_name")) {
        cout << "algorithm_name required" << endl;

        cout << "see --help for more information" << endl;
    } else{
        cout << vm["algorithm_name"].as<string>() << endl;
        setvm(vm);
        addArg("dataset", (const float**)data);
        addArg("queryset", (const float**)query);
        addArg("ground_truth", (const Result**)results);
        
        run(vm["algorithm_name"].as<string>());
    }
    return 0;
}