// ***************************************************************************
// ParaModel.h (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _PARAMODEL_H
#define _PARAMODEL_H

#include "bayesopt/bayesopt.hpp"
#include "specialtypes.hpp"

class ModelTwoD: public bayesopt::ContinuousModel {
	public:
		ModelTwoD(bayesopt::Parameters paras):
			ContinuousModel(2, paras) {}

		double evaluateSample(const vectord& paras);
		bool checkReachability(const vectord &query);
		
};

class ParaModel {
	public:
		void init(bayesopt::Parameters &paras);
		void learn();
};

#endif
