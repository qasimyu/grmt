// ***************************************************************************
// InputParser.cpp (c) 2020 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <limits.h>

#include "InputParser.h"
#include "MyDefine.h"

using namespace std;

void InputParser::parseArgs(int argc, char *argv[]) {
	string inputFile = "", rinputFile = "", outputPrefix = "";
	string clabelFile = "", mlabelFile = "";
	int maxl = 0, beamsize = 5, threads = 1;
	double lambda = 0.7, kappa = 1;
	double beta = -1, alpha = -1;
	int n_init = 100, n_iter = 30;
	
	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"input", required_argument, 0, 'i'},
		//{"rinput", required_argument, 0, 'r'},
		{"output", required_argument, 0, 'o'},
		{"clabel", required_argument, 0, 'c'},
		{"mlabel", required_argument, 0, 'm'},
		{"lambda", required_argument, 0, 'l'},
		{"kappa", required_argument, 0, 'K'},
		//{"bsize", required_argument, 0, 'n'},
		{"maxl", required_argument, 0, 'k'},
		{"threads", required_argument, 0, 't'},
		{"alpha", required_argument, 0, 'a'},
		{"beta", required_argument, 0, 'b'},
		{"n_init", required_argument, 0, 'n'},
		{"n_iter", required_argument, 0, 'N'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hvi:o:c:m:l:K:k:t:a:b:n:N:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage(argv[0]);
				exit(0);
			case 'v':
				cerr << "GRMT version " << current_version << endl;
				exit(0);
			case 'i':
				inputFile = optarg;
				break;
			/*
			case 'r':
				rinputFile = optarg;
				break;
			*/
			case 'o':
				outputPrefix = optarg;
				break;
			case 'c':
				clabelFile = optarg;
				break;
			case 'm':
				mlabelFile = optarg;
				break;
			case 'l':
				lambda = atof(optarg);
				break;
			case 'K':
				kappa = atof(optarg);
				break;
			/*
			case 'n':
				beamsize = atoi(optarg);
				break;
			*/
			case 'k':
				maxl = atoi(optarg);
				break;
			case 't':
				threads = atoi(optarg);
				break;
			case 'a':
				alpha = atof(optarg);
				break;
			case 'b':
				beta = atof(optarg);
				break;
			case 'n':
				n_init = atoi(optarg);
				break;
			case 'N':
				n_iter = atoi(optarg);
				break;
			default :
				usage(argv[0]);
				exit(1);
		}
	}
	
	if(inputFile.empty()){
        cerr << "Use --input to specify the file containing mutation data." << endl;
		usage(argv[0]);
        exit(1);
    }
	/*
	if(rinputFile.empty()){
        cerr << "Warning: the file containing real mutation data was not specified." << endl;
    }
	*/

	if(outputPrefix.empty()){
		cerr << "Use --output to specify the prefix of result file names." << endl;
		usage(argv[0]);
		exit(1);
	}
	/*
	if(beamsize < 1) {
		cerr << "Error: the value of parameter \"bsize\" should be a positive integer." << endl;
		usage(argv[0]);
		exit(1);
	}
	*/
	
	if(lambda < 0.5 || lambda > 1) {
		cerr << "Error: the value of parameter \"lambda\" should be in (0.5, 1]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(kappa <= 0) {
		cerr << "Error: the value of parameter \"kappa\" should be a positive number." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(maxl < 0) {
		cerr << "Error: the value of parameter \"maxl\" should be a non-negative integer." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(threads < 1) {
		cerr << "Error: the value of parameter \"threads\" should be a positive integer." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(alpha > 1) {
		cerr << "Error: the value of parameter \"alpha\" should be in [0, 1]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(beta > 1) {
		cerr << "Error: the value of parameter \"beta\" should be in [0, 1]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(n_init < 10) {
		cerr << "Error: the value of parameter \"n_init\" should be in [10,]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(n_iter < 10) {
		cerr << "Error: the value of parameter \"n_iter\" should be in [10,]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	config.setStringPara("input", inputFile);
	//config.setStringPara("rinput", rinputFile);
	config.setStringPara("output", outputPrefix);
	config.setStringPara("clabel", clabelFile);
	config.setStringPara("mlabel", mlabelFile);
	//config.setIntPara("beamsize", beamsize);
	config.setIntPara("maxl", maxl);
	config.setIntPara("threads", threads);
	config.setIntPara("n_init", n_init);
	config.setIntPara("n_iter", n_iter);
	config.setRealPara("lambda", lambda);
	config.setRealPara("kappa", kappa);
	config.setRealPara("alpha", alpha);
	config.setRealPara("beta", beta);
	
	/*** check output directory ***/
	size_t i = outputPrefix.find_last_of('/');
	string outputDir = outputPrefix;
	if(i != string::npos) {
		outputDir = outputPrefix.substr(0, i);
	}
	bool is_exist = (access(outputDir.c_str(), F_OK) == 0);
	if(!is_exist) {
		cerr << "Error: the output directory " << outputDir << " does not exist!" << endl;
		exit(1);
	}
	
	/*** create thread pool ***/
	threadpool = new ThreadPool(threads);
	threadpool->pool_init();
	
	string binaryPath = argv[0];
	char abs_path_buff[PATH_MAX];
	realpath(binaryPath.c_str(), abs_path_buff);
	string binPath = abs_path_buff;
	int indx = binPath.rfind("bin");
	binPath = binPath.substr(0, indx-1);
	config.setStringPara("binPath", binPath);
	
}

void InputParser::usage(const char* app) {
	cerr << "Usage: " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -v, --version                   print software version" << endl
		<< "    -i, --input <string>            mutation data file" << endl
		//<< "    -r, --rinput <string>           input file containing real mutation data" << endl
		<< "    -o, --output <string>           prefix of output file names" << endl
		<< "    -c, --clabel <string>           file defining labels of the cells" << endl
		<< "    -m, --mlabel <string>           file defining labels of the mutations" << endl
		<< "    -l, --lambda <double>           value of hyper-parameter lambda [default:0.7]" << endl
		<< "    -K, --kappa <double>            value of hyper-parameter kappa [default:1]" << endl
		//<< "    -n, --bsize <int>               beam search size [default:5]" << endl
		<< "    -k, --maxl <int>                maximum number of times that a mutation can be lost [default:0]" << endl
		<< "    -t, --threads <int>             number of threads to use [default:1]" << endl
		<< "    -a, --alpha <double>            false positive rate" << endl
		<< "    -b, --beta <double>             false negative rate" << endl
		<< "    -n, --n_init <int>              number of initial points to sample in the BO algorithm [default:100]" << endl
		<< "    -N, --n_iter <int>              number of points to evaluate in the BO algorithm [default:30]" << endl
		<< endl
		<< "Example:" << endl
		<< app << " -i ./testdata/example.txt -o ./testdata/example -a 0.01 -b 0.2 " << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}
