// ***************************************************************************
// ParaModel.cpp (c) 2020 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cmath>

#include "ParaModel.h"
#include "Matrix.h"
#include "MyDefine.h"

using namespace std;

double ModelTwoD::evaluateSample(const vectord& paras) {
	if(paras.size() < 2) {
		cerr << "Error: no enough inputs for 2D model." << endl;
		return 0;
	}
	if(paras.size() != 2) {
		cerr << "WARNING: This only works for 2D inputs." << endl
			<< "Using only first two components." << endl;
	}
	return -treecaller.evaluate(paras(0), paras(1));
}

bool ModelTwoD::checkReachability(const vectord &query) {
	return true;
}

void ParaModel::init(bayesopt::Parameters &paras) {
	paras.n_init_samples = config.getIntPara("n_init");
	paras.n_iterations = config.getIntPara("n_iter");
	paras.verbose_level = 0; // 1 for info
	paras.random_seed = 0;
//	paras.init_method = 0;
}

void ParaModel::learn() {	
	double alpha = config.getRealPara("alpha");
	double beta = config.getRealPara("beta");
	double max_beta = config.getRealPara("max_beta");
	
	if(alpha < 0 || beta < 0) {
		bayesopt::Parameters paras;
		init(paras);
		ModelTwoD model(paras);
		
		vectord paraLowerBound(2), paraUpperBound(2);
		if(alpha >= 0) {
			paraLowerBound(0) = alpha;
			paraUpperBound(0) = alpha;
		}
		else {
			paraLowerBound(0) = 5e-6;//5e-6
			paraUpperBound(0) = 0.05;//0.05
		}
		if(beta >= 0) {
			paraLowerBound(1) = beta;
			paraUpperBound(1) = beta;
		}
		else {
			paraLowerBound(1) = 0.01;
			paraUpperBound(1) = max_beta;
		}
		model.setBoundingBox(paraLowerBound, paraUpperBound);
		
		vectord result(2);
		model.optimize(result);
		treecaller.setParas(result(0), result(1));
	}
	else {
		treecaller.evaluate(alpha, beta);
	}
	
	treecaller.saveResults();
	
}
