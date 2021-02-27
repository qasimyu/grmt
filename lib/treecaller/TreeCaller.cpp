// ***************************************************************************
// TreeCaller.cpp (c) 2020 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <limits>

#include "TreeCaller.h"
#include "split.h"
#include "MyDefine.h"

using namespace std;


TreeCaller::TreeCaller() {
}

void TreeCaller::loadData() {
	loadMutationData();
	//loadRealData();
	loadCellLabels();
	loadMutaLabels();
}

void TreeCaller::loadMutationData() {
	string inputFile = config.getStringPara("input");
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	fgets(buf, 500, fp);
	fclose(fp);
	int num_cell = atoi(buf);
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_muta = 0;
	missing_rate = 0;
	homo_muta = false;
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_muta == 0) {
			num_muta = fields.size();
			obsData.resize(num_cell, num_muta, false);
		}
		if(num_muta != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_muta; i++) {
			indx = (line_num-1)*num_muta+i;
			k = atoi(fields[i].c_str());
			if(k == 3) {
				missing_rate++;
			}
			if(k == 2) {
				homo_muta = true;
			}
			obsData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_cell);
	missing_rate /= (num_cell*num_muta);
	
	cerr << "Total " << num_cell << " cells and " << num_muta << " mutations were loaded from file " << inputFile << endl;
}

void TreeCaller::loadRealData() {
	string inputFile = config.getStringPara("rinput");
	if(inputFile.empty()) {
		return;
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	fgets(buf, 500, fp);
	fclose(fp);
	int num_cell = atoi(buf)-1;
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_muta = 0;
	
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	getline(ifs, line);
	vector<string> fields = split(line, ':');
	if(fields.size() > 1) {
		fields = split(fields[1], '\t');
		for(i = 0; i < fields.size(); i++) {
			doublet_indxs.push_back(atoi(fields[i].c_str())-1);
		}
	}
	
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_muta == 0) {
			num_muta = fields.size();
			realData.resize(num_cell, num_muta, false);
		}
		if(num_muta != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_muta; i++) {
			indx = (line_num-1)*num_muta+i;
			k = atoi(fields[i].c_str());
			realData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_cell);
	
	cerr << "real mutation data were loaded from file " << inputFile << endl;
}

void TreeCaller::loadCellLabels() {
	string inputFile = config.getStringPara("clabel");
	if(inputFile.empty()) {
		char buf[100];
		for(int i = 0; i < obsData.getROWS(); i++) {
			sprintf(buf, "%d", i+1);
			cLabels.push_back(buf);
		}
		return;
	}
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	long line_num = 0;
	string line;
	while(getline(ifs, line)) {
		line_num++;
		if(line.empty()) {
			cerr << "Error: malformed label file " << inputFile <<
                    ", empty item is not allowed @line " << line_num << endl;
			exit(1);
		}
		cLabels.push_back(line);
	}
	ifs.close();
	
	if(cLabels.size() != obsData.getROWS()) {
		cerr << "Error: the number of cell labels is not consistent with the number of cells." << endl;
		exit(1);
	}
	
	cerr << "the labels of cells were loaded from file " << inputFile << endl;
}

void TreeCaller::loadMutaLabels() {
	string inputFile = config.getStringPara("mlabel");
	if(inputFile.empty()) {
		char buf[100];
		for(int i = 1; i <= obsData.getCOLS(); i++) {
			sprintf(buf, "%d", i);
			mLabels.push_back(buf);
		}
		return;
	}
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	long line_num = 0;
	string line;
	while(getline(ifs, line)) {
		line_num++;
		if(line.empty()) {
			cerr << "Error: malformed label file " << inputFile <<
                    ", empty item is not allowed @line " << line_num << endl;
			exit(1);
		}
		mLabels.push_back(line);
	}
	ifs.close();
	
	if(mLabels.size() != obsData.getCOLS()) {
		cerr << "Error: the number of mutation labels is not consistent with the number of mutations." << endl;
		exit(1);
	}
	
	cerr << "the labels of mutations were loaded from file " << inputFile << endl;
}

void TreeCaller::initParas() {
	int i, j;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	int k = config.getIntPara("maxl"); // k-Dollo evolutionary model
	int num_node = 1+num_muta*(k+1);
	
	Tree tree(-1, -1, 0);
	Matrix<int> trees(1, num_node, -1); //initial trees, only contains root node
	tree.trees = trees;
	
	Matrix<int> nodeStates(num_node, num_muta, -1);
	nodeStates.setRow(0, 0);
	tree.nodeStates = nodeStates;
	Matrix<int> cellLocs(1, num_cell, 0);
	tree.cellLocs = cellLocs;
	Matrix<int> pathLens(1, num_node, 0);
	tree.pathLens = pathLens;
	
	double eps = numeric_limits<double>::epsilon();
	double score = 0;
	Matrix<double> obsLLs(num_cell, num_node+1, (double) numeric_limits<long>::min());
	for(i = 0; i < num_cell; i++) {
		double ll = 0;
		for(j = 0; j < num_muta; j++) {
			ll += getLL(0, obsData[i*num_muta+j]);
		}
		score += ll;
		obsLLs[i*(num_node+1)] = ll;
		obsLLs[i*(num_node+1)+num_node] = ll;
	}
	tree.score = score;
	tree.obsLLs = obsLLs;
	tree.singlet_count = 1;
	
	bestTree = tree;
	bs_Trees.clear();
	bs_Trees.push_back(tree);
}

double TreeCaller::evaluate(double alpha, double beta) {
	setParas(alpha, beta);
	call();
	return bestTree.score;
}

void TreeCaller::call() {
	int i, j;
	
	//initialize trees and parameters
	initParas();
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	
	int k = config.getIntPara("maxl"); // k-Dollo evolutionary model
	int num_node = 1+num_muta*(k+1);
	
	double score_pre, acc;
	
	int iter = 0;
	while(1) {
		score_pre = bestTree.score;
		
		findAndEvaluateNeighbors();
		
		if(candi_trees_all.empty()) { //search finished
			break;
		}
		
		//find top ranked trees
		findTopRankedTrees();
		
		if(bs_Trees.empty()) {
			break;
		}
		
		iter++;
		
		//if(iter%100 == 0) {
			
			//cerr << "--------------- search report -----------------" << endl;
			//cerr << "iteration " << iter << " done." << endl;
			//cerr << "best score: " << bestTree.score << endl;
			/*
			acc = evalAcc(bestTree);
			if(acc == -1) {
				cerr << "accuracy: NA" << endl;
			}
			else {
				cerr << "accuracy: " << acc << endl;
			}
			*/
			//cerr << "best tree: ";
			//printTree(bestTree);
			/*
			cerr << "obsProbs:" << endl;
			for(i = 0; i < bs_Trees.size(); i++) {
				bs_Trees[i].obsLLs.Print();
			}
			//cerr << endl;
			*/
			
			//cerr << "cellLocs: ";
			//bestTree.cellLocs.Print();
			/*
			cerr << "cellLLs: " << endl;
			bestTree.obsLLs.Print();
			
			
			cerr << "scores:";
			for(i = 0; i < bs_Trees.size(); i++) {
				cerr << " " << bs_Trees[i].score;
			}
			cerr << endl;
			
			cerr << "weight:";
			for(i = 0; i < bs_Trees.size(); i++) {
				cerr << " " << bs_Trees[i].weight;
			}
			cerr << endl;
			
			cerr << "accuracies:";
			for(i = 0; i < bs_Trees.size(); i++) {
				cerr << " " << evalAcc(bs_Trees[i]);
			}
			cerr << endl;

			cerr << "trees:" << endl;
			for(i = 0; i < bs_Trees.size(); i++) {
				printTree(bs_Trees[i]);
			}
			*/		
			//cerr << "--------------- search report -----------------" << endl;
		//}
	}
	
	cerr << "--------------- search report -----------------" << endl;
	cerr << "alpha: " << alpha << endl;
	cerr << "beta: " << beta << endl;
	cerr << "score: " << bestTree.score << endl;
	
	/*
	acc = evalAcc(bestTree);
	if(acc == -1) {
		cerr << "accuracy: NA" << endl;
	}
	else {
		cerr << "accuracy: " << acc << endl;
	}
	*/
	//cerr << "tree: ";
	//printTree(bestTree);
	cerr << "--------------- search report -----------------" << endl;
	
	Solution s(alpha, beta, bestTree);
	solutions.push_back(s);
}

double TreeCaller::evalAcc(Tree& tree) {
	if(realData.getROWS() == 0) {
		return -1;
	}
	int i, j, k;
	int r, c;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	int num_node = bestTree.trees.getCOLS();
	
	//Matrix<int> cellLocs;
	//inferCellLocs(tree, cellLocs);
	
	Matrix<int>& cellLocs = tree.cellLocs;
	Matrix<int>& nodeStates = tree.nodeStates;
	
	double acc = 0;
	int n = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < doublet_indxs.size(); j++) {
			if(doublet_indxs[j] == i) {
				break;
			}
		}
		if(j < doublet_indxs.size()) {
			continue;
		}
		n++;
		
		k = cellLocs[i];
		for(j = 0; j < num_muta; j++) {
			if(nodeStates[k*num_muta+j] == realData[i*num_muta+j]) {
				acc++;
			}
		}
	}
	acc /= (n*num_muta);
	
	return acc;
}

void TreeCaller::printTree(Tree& tree) {
	tree.trees.Print();
}

void TreeCaller::findTopRankedTrees() {
	int i, j, k, m;
	//int bs = config.getIntPara("beamsize");
	int bs = 1;
	int num_tree = bs_Trees.size();
	int mindx;
	
	if(bs_Trees.size() > 1) {
		double max_weight = 0;
		vector<int> indxs;
		for(i = 0; i < bs_Trees.size(); i++) {
			vector<Tree>& candi_trees = candi_trees_all[i];
			double weight, weight_m = 0;
			int mindx1;
			for(m = 0; m < candi_trees.size(); m++) {
				weight = candi_trees[m].weight;
				if(weight_m < weight) {
					weight_m = weight;
					mindx1 = m;
				}
			}
			if(max_weight < weight_m) {
				max_weight = weight_m;
				mindx = i;
			}
			updateTree(candi_trees[mindx1]);
			indxs.push_back(mindx1);
		}
		
		bs_Trees.clear();
		for(i = 0; i < indxs.size(); i++) {
			j = indxs[i];
			bs_Trees.push_back(candi_trees_all[i][j]);
		}
	}
	else {
		//map<int, vector<int> > sel_indxs;
		//for(i = 0; i < num_tree; i++) {
			vector<Tree>& candi_trees = candi_trees_all[0];
			Matrix<double> topWeights(1, bs, (double) numeric_limits<long>::min());
			Matrix<int> topIndxs(1, bs, -1);
			topWeights[0] = candi_trees[0].weight;
			topIndxs[0] = 0;
			for(m = 1; m < candi_trees.size(); m++) {
				double weight = candi_trees[m].weight;
				for(j = 0; j < bs; j++) {
					if(weight > topWeights[j]) {
						for(k = bs-1; k > j; k--) {
							topWeights[k] = topWeights[k-1];
							topIndxs[k] = topIndxs[k-1];
						}
						topWeights[j] = weight;
						topIndxs[j] = m;
						break;
					}
				}
			}
			
			for(j = 0; j < bs; j++) {
				k = topIndxs[j];
				if(k != -1) {
					updateTree(candi_trees[k]);
				}
			}
			
		//}
		
		bs_Trees.clear();
		for(j = 0; j < bs; j++) {
			k = topIndxs[j];
			if(k != -1) {
				bs_Trees.push_back(candi_trees[k]);
			}
		}
		mindx = 0;
		
	}
	
	double kappa = config.getRealPara("kappa");
	
	if(bs_Trees[mindx].weight < kappa) {
		bs_Trees.clear();
		return;
	}
	
	if(bs_Trees[mindx].score > bestTree.score) {
		bestTree = bs_Trees[mindx];
		//cerr << bs_Trees[mindx].score << endl;
		//cerr << bestTree.parNode << "->" << bestTree.childNode << endl;
	}
}

void TreeCaller::updateTree(Tree& tree) {
	int i, j, c, n;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	int num_node = bs_Trees[0].trees.getCOLS();
	int k = config.getIntPara("maxl");
	int parTree = tree.parTree;
	int parNode = tree.parNode;
	int childNode = tree.childNode;
	
	tree.trees = bs_Trees[parTree].trees;
	tree.trees[childNode] = parNode;
	tree.pathLens = bs_Trees[parTree].pathLens;
	tree.pathLens[childNode] = tree.pathLens[parNode]+1;
	Matrix<int>& pathLens = tree.pathLens;
	
	Matrix<int> nodeStates = bs_Trees[parTree].nodeStates;
	Matrix<int> cellLocs = bs_Trees[parTree].cellLocs;
	Matrix<double> obsLLs = bs_Trees[parTree].obsLLs;
	
	int mut_indx = (childNode-1)/(k+1); // the index of mutation
	int m = (childNode-1)%(k+1); // m=0 represents mutation gain, other values denotes mutation loss
	
	for(j = 0; j < num_muta; j++) {
		nodeStates[childNode*num_muta+j] = nodeStates[parNode*num_muta+j];
	}
	if(m == 0) { // gain mutation
		nodeStates[childNode*num_muta+mut_indx] = 1;
	}
	else { // lose mutation
		nodeStates[childNode*num_muta+mut_indx] = 0;
	}
	tree.nodeStates = nodeStates;
	
	Matrix<int>& trees = bs_Trees[parTree].trees;
	for(c = 0; c < num_cell; c++) {
		j = cellLocs[c];
		double ll_p = obsLLs[c*(num_node+1)+parNode];
		double tmp1 = getLL(1-nodeStates[childNode*num_muta+mut_indx], obsData[c*num_muta+mut_indx]);
		double tmp2 = getLL(nodeStates[childNode*num_muta+mut_indx], obsData[c*num_muta+mut_indx]);
		double ll_c = ll_p-tmp1+tmp2;
		obsLLs[c*(num_node+1)+childNode] = ll_c;
		obsLLs[c*(num_node+1)+num_node] = log(exp(obsLLs[c*(num_node+1)+num_node])+exp(ll_c));
		
		if(obsLLs[c*(num_node+1)+j] < ll_c) {
			cellLocs[c] = childNode;
		}
	}
	tree.cellLocs = cellLocs;
	tree.obsLLs = obsLLs;
	
}

void TreeCaller::findAndEvaluateNeighbors() {
	int i, j, n, ii;
	int num_node = bs_Trees[0].trees.getCOLS();
	int k = config.getIntPara("maxl");
	int num_muta = obsData.getCOLS();
	
	candi_trees_all.clear();
	
	for(n = 0; n < bs_Trees.size(); n++) {
		Matrix<int>& trees = bs_Trees[n].trees;
		Matrix<int>& nodeStates = bs_Trees[n].nodeStates;
		bool init_tumor = false;
		for(i = 1; i < num_node; i++) {
			if(trees[i] == 0) {
				init_tumor = true;
				break;
			}
		}
		for(i = 1; i < num_node; i++) {
			if(trees[i] != -1) { // node i is already added to the tree
				continue;
			}
			int mut_indx = (i-1)/(k+1); // the index of mutation
			int m = (i-1)%(k+1); // m=0 represents mutation gain, other values denotes mutation loss
			//attach node i to one node of the tree
			for(j = 0; j < num_node; j++) {
				if(j != 0 && trees[j] == -1) { // node j is not yet in the tree
					continue;
				}
				if(j == 0 && init_tumor) {
					continue;
				}
				if(j > 0 && (j-1)/(k+1) == mut_indx) {
					continue;
				}
				if(m != 0 && nodeStates[j*num_muta+mut_indx] == 0) {
					continue;
				}
				
				Tree tree(n, j, i);
				candi_trees_all[n].push_back(tree);
			}
		}
	}
	if(!candi_trees_all.empty()) {
		int threads = threadpool->getThreadNumber();
		vector<int*> threadParas;
		for(n = 0; n < bs_Trees.size(); n++) {
			int N = max(100, (int) candi_trees_all[n].size()/threads);
			i = 0;
			while(i < candi_trees_all[n].size()) {
				int *tparas = new int[3];
				tparas[0] = n; tparas[1] = i;
				if(i+N <= candi_trees_all[n].size()) {
					tparas[2] = i+N-1;
				}
				else {
					tparas[2] = candi_trees_all[n].size()-1;
				}
				threadpool->pool_add_work(&TreeCaller::evaluateCandiTrees, tparas, 0);
				threadParas.push_back(tparas);
				i += N;
			}
		}
		threadpool->wait();
		for(i = 0; i < threadParas.size(); i++) {
			delete[] threadParas[i];
		}
	}
	
}

void* TreeCaller::evaluateCandiTrees(const void *arg) {
	int *tmp = (int*) arg;
	int tindx = tmp[0];
	int sindx = tmp[1];
	int eindx = tmp[2];
	int i, j, c, n, m;
	int mut_indx;
	
	double eps = numeric_limits<double>::epsilon();
	
	Matrix<int>& obsData = treecaller.getObsData();
	
	vector<Tree>& candi_trees_all = treecaller.getCandiTrees(tindx);
	vector<Tree>& bs_Trees = treecaller.getTrees();
	
	int num_node = bs_Trees[0].trees.getCOLS();
	
	double lambda = config.getRealPara("lambda");
	int k = config.getIntPara("maxl");
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	Matrix<int> states(1, num_muta);
	
	for(i = sindx; i <= eindx; i++) {
		Tree& tree = candi_trees_all[i];
		int parTree = tree.parTree;
		int parNode = tree.parNode;
		int childNode = tree.childNode;
		
		Tree& par_tree = bs_Trees[parTree];
		
		Matrix<int>& trees = par_tree.trees;
		Matrix<int>& nodeStates = par_tree.nodeStates;
		Matrix<int>& cellLocs = par_tree.cellLocs;
		Matrix<int>& pathLens = par_tree.pathLens;
		Matrix<double>& obsLLs = par_tree.obsLLs;
		
		mut_indx = (childNode-1)/(k+1); // the index of mutation
		m = (childNode-1)%(k+1); // m=0 represents mutation gain, other values denotes mutation loss
		
		for(j = 0; j < num_muta; j++) {
			states[j] = nodeStates[parNode*num_muta+j];
		}
		if(m == 0) { // gain mutation
			states[mut_indx] = 1;
		}
		else { // lose mutation
			states[mut_indx] = 0;
		}
		double count_pre = 0, count_p = 0, count = 0, score = 0;
		int singlet_count = par_tree.singlet_count+1;
		for(c = 0; c < num_cell; c++) {
			j = cellLocs[c];
			double ll_p = obsLLs[c*(num_node+1)+parNode];
			count_pre += exp(ll_p-obsLLs[c*(num_node+1)+num_node]);
			double tmp1 = treecaller.getLL(1-states[mut_indx], obsData[c*num_muta+mut_indx]);
			double tmp2 = treecaller.getLL(states[mut_indx], obsData[c*num_muta+mut_indx]);
			double ll_c = ll_p-tmp1+tmp2;
			double prob_singlet = exp(obsLLs[c*(num_node+1)+num_node])+exp(ll_c);
			
			count += exp(ll_c-log(prob_singlet));
			count_p += exp(ll_p-log(prob_singlet));
			score += max(obsLLs[c*(num_node+1)+j], ll_c);
			//score += log(prob_singlet/singlet_count);
		}
		tree.weight = (count-count_pre+count_p)*(1-lambda)+(count_pre-count_p)*lambda;
		tree.score = score;
		tree.cell_attached = count;
		tree.singlet_count = singlet_count;
		//cerr << parNode << ", " << childNode << ", " << count_pre-count_p << ", " << count << ", " << tree.weight << endl;
	}
	return NULL;
}

void TreeCaller::inferCellLocs(Tree& tree, Matrix<int>& cellLocs) {
	int i, j, k;
	int num_node = tree.trees.getCOLS();
	int num_cell = obsData.getROWS();
	Matrix<double>& obsLLs = tree.obsLLs;
	Matrix<int>& pathLens = tree.pathLens;
	cellLocs.resize(1, num_cell, false);
	
	for(i = 0; i < num_cell; i++) {
		double max_ll = obsLLs[i*(num_node+1)];
		for(j = 1; j < num_node; j++) {
			if(obsLLs[i*(num_node+1)+j] > max_ll) {
				max_ll = obsLLs[i*(num_node+1)+j];
			}
		}
		int min_pl = numeric_limits<int>::max();
		for(j = 0; j < num_node; j++) {
			if(abs(obsLLs[i*(num_node+1)+j]-max_ll) < 1.0e-5) {
				if(pathLens[j] < min_pl) {
					min_pl = pathLens[j];
					k = j;
				}
			}
		}
		cellLocs[i] = k;
	}
	
}

double TreeCaller::getLL(int gt, int observed) {
	double prob = getObsProb(gt, observed);
	return log(prob);
}

double TreeCaller::getObsProb(int gt, int observed) {
	double prob = 1, eps = numeric_limits<double>::epsilon();
	if(gt == 1) {
		if(observed == 0) {
			prob = (homo_muta)? beta/2+eps:beta+eps;
		}
		else if(observed == 1) {
			prob = 1-beta+eps;
		}
		else if(observed == 2) {
			prob = beta/2+eps;
		}
		/*
		else {
			prob = missing_rate/(1-missing_rate);
		}
		prob *= (1-missing_rate);
		*/
	}
	else {
		if(observed == 0) {
			prob = (homo_muta)? 1-alpha-alpha*beta/2+eps:1-alpha+eps;
		}
		else if(observed == 1) {
			prob = alpha+eps;
		}
		else if(observed == 2) {
			prob = alpha*beta/2+eps;
		}
		/*
		else {
			prob = missing_rate/(1-missing_rate);
		}
		prob *= (1-missing_rate);
		*/
		
	}
	return prob;
}

void TreeCaller::saveResults() {
	int i, j, m;
	int k = config.getIntPara("maxl");
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	
	int best_s_indx = 0;
	double maxScore = solutions[0].tree.score;
	for(i = 1; i < solutions.size(); i++) {
		if(solutions[i].tree.score > maxScore) {
			best_s_indx = i;
			maxScore = solutions[i].tree.score;
		}
	}
	alpha = solutions[best_s_indx].alpha;
	beta = solutions[best_s_indx].beta;
	bestTree = solutions[best_s_indx].tree;
	
	cerr << "--------------- best solution -----------------" << endl;
	cerr << "alpha: " << alpha << endl;
	cerr << "beta: " << beta << endl;
	cerr << "score: " << bestTree.score << endl;
	/*
	double acc = evalAcc(bestTree);
	if(acc == -1) {
		cerr << "accuracy: NA" << endl;
	}
	else {
		cerr << "accuracy: " << acc << endl;
	}
	*/
	cerr << "-----------------------------------------------" << endl;
	
	//Matrix<int> cellLocs = bestTree.cellLocs;
	Matrix<int> cellLocs;
	inferCellLocs(bestTree, cellLocs);
	//cellLocs.Print();
	
	Matrix<int>& trees = bestTree.trees;
	int num_node = trees.getCOLS();
	
	Matrix<int> cell_attached(1, num_node, 0);
	Matrix<int> child_counts(1, num_node, 0);
	for(i = 0; i < num_cell; i++) {
		j = cellLocs[i];
		cell_attached[j] = 1;
	}
	for(i = 1; i < num_node; i++) {
		j = trees[i];
		if(j >= 0) {
			child_counts[j] += 1;
		}
	}
	for(i = 0; i < num_node; i++) {
		if(child_counts[i] == 0 && cell_attached[i] == 0) {
			trees[i] = -1;
		}
	}
	
	string outputPrefix = config.getStringPara("output");
	
	string fn = outputPrefix+".solution.txt";
	ofstream ofs;
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	ofs << "ll = " << bestTree.score << endl;
	ofs << "alpha = " << alpha << endl;
	ofs << "beta = " << beta << endl;
	/*
	ofs << "tree = ";
	for(i = 0; i < num_node; i++) {
		if(i < num_node-1)
			ofs << trees[i] << " ";
		else
			ofs << trees[i] << endl;
	}
	*/
	ofs.close();
	
	fn = outputPrefix+".dot";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	
	int mut_indx;
	ofs << "digraph T {" << endl;
	ofs << "node [color=deeppink4, style=filled, fontcolor=white];" << endl;
	ofs << "0 [label=\"\"];" << endl;
	for(i = 1; i < num_node; i++) {
		if(trees[i] == -1) {
			continue;
		}
		mut_indx = (i-1)/(k+1);
		m = (i-1)%(k+1);
		if(m == 0) {
			ofs << i << " [label=\"" << mLabels[mut_indx] << "\"];" << endl;
		}
		else {
			ofs << i << " [color=\".7 .3 1.0\",label=\"-" << mLabels[mut_indx] << "\"];" << endl;
		}
	}
	
	map<int, string> ulabels;
	map<int, int> counts;
	for(i = 0; i < num_cell; i++) {
		int key = cellLocs[i];
		if(ulabels[key].empty()) {
			ulabels[key] = cLabels[i];
			counts[key] = 0;
		}
		else {
			counts[key]++;
			if(counts[key] == 5) {
				counts[key] = 0;
				ulabels[key] = ulabels[key]+" "+cLabels[i]+"\\n";
			}
			else {
				ulabels[key] = ulabels[key]+" "+cLabels[i];
			}
		}
	}
	ofs << "node [color=lightgrey, style=filled, fontcolor=black];" << endl;
	map<int, string>::iterator it;
	for(i = num_node, it = ulabels.begin(); it != ulabels.end(); i++, it++) {
		ofs << i << " [label=\"" << it->second << "\"];" << endl;
	}
	
	for(i = 0; i < num_node; i++) {
		if(trees[i] == -1) {
			continue;
		}	
		ofs << trees[i] << "->" << i << ";" << endl;
	}
	for(i = num_node, it = ulabels.begin(); it != ulabels.end(); i++, it++) {
		ofs << it->first << "->" << i << ";" << endl;
	}
	ofs << "}" << endl;
	
	ofs.close();
	
}

