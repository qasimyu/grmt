// ***************************************************************************
// TreeCaller.h (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _TREECALLER_H
#define _TREECALLER_H

#include <map>
#include <vector>
#include <pthread.h>

#include "Matrix.h"

class Tree {
	public:
		Tree() {}
		Tree(int parTree, int parNode, int childNode) :
			parTree(parTree), parNode(parNode), childNode(childNode) {}
		int parTree;
		int parNode;
		int childNode;
		double score, weight;
		double cell_attached;
		int singlet_count;
		
		Matrix<int> trees;
		Matrix<int> nodeStates, cellLocs;
		Matrix<int> pathLens;
		Matrix<double> obsLLs;
};

class Solution {
	public:
		Solution(double alpha, double beta, Tree& tree) :
			alpha(alpha), beta(beta), tree(tree) {}
		double alpha;
		double beta;
		Tree tree;
};

class TreeCaller {
	private:
		Matrix<int> obsData;
		Matrix<int> realData;
		vector<string> cLabels;
		vector<string> mLabels;
		bool homo_muta;
		double alpha, beta;
		
		double missing_rate;
		
		vector<int> doublet_indxs;
		
		vector<Solution> solutions;
		map<int, vector<Tree> > candi_trees_all;
		
		vector<Tree> bs_Trees;
		
		Tree bestTree;
		
		void loadMutationData();
		void loadRealData();
		void loadCellLabels();
		void loadMutaLabels();
		
		void initParas();
		void findTopRankedTrees();
		static void* evaluateCandiTrees(const void *arg);
		
		void printTree(Tree& tree);
		void updateTree(Tree& tree);
		
		void findAndEvaluateNeighbors();
		void inferCellLocs(Tree& tree, Matrix<int>& cellLocs);
		
		double getLL(int gt, int observed);
		double getObsProb(int gt, int observed);
		
		double evalAcc(Tree& tree);
	public:
		TreeCaller();
		
		void loadData();
		Matrix<int>& getObsData() {return obsData;}
		
		vector<Tree>& getCandiTrees(int tindx) {return candi_trees_all[tindx];}
		
		vector<Tree>& getTrees() {return bs_Trees;}
		
		void setParas(double alpha, double beta) {this->alpha = alpha; this->beta = beta;}
		
		double evaluate(double alpha, double beta);
		
		void call();
		void saveResults();
};

#endif
