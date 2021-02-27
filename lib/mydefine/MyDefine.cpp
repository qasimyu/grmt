// ***************************************************************************
// MyDefine.cpp (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdio>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

#include "MyDefine.h"

using namespace std;

MyDefine::MyDefine() {
}

MyDefine::~MyDefine() {
}

string current_version = "1.0";

/*** definition of global vars ***/
Config config;
InputParser parser;
TreeCaller treecaller;
ParaModel paramodel;
ThreadPool *threadpool;
/*** end of definition ***/

double normrand(double mu, double sigma) {
	double eps = numeric_limits<double>::epsilon();
    double pi = 3.14159265358979323846;

    static double z0, z1;
    static bool flag = true;

    if(!flag) {
       return z1*sigma+mu;
	}
	flag = !flag;
	
    double u1, u2;
    do {
       u1 = rand()*(1.0/RAND_MAX);
       u2 = rand()*(1.0/RAND_MAX);
	}while(u1 <= eps);

    z0 = sqrt(-2.0*log(u1))*cos(2*pi*u2);
    z1 = sqrt(-2.0*log(u1))*sin(2*pi*u2);
	
    return z0*sigma+mu;
}

//****** trim string ******//
string trim(const string &str, const char *charlist) {
	string ret(str);
	size_t indx = ret.find_first_not_of(charlist);
	if(indx != string::npos) {
		ret.erase(0, indx);
		indx = ret.find_last_not_of(charlist);
		ret.erase(indx+1);
	}
	else {
		ret.erase();
	}
	return ret;
}

//****** produce random double ******//
double randomDouble(double start, double end) {
	return start+(end-start)*(rand()/(RAND_MAX+1.0));
}

//****** produce random int ******//
long randomInteger(long start, long end) {
	return start+(end-start)*(rand()/(RAND_MAX+1.0));
}



