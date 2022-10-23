#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <vector>
#include <string>
#include <float.h>
#include <omp.h>
#include "utils.h"

using namespace std;
using namespace arma;

double RMI(arma::Col<int> &g1, arma::Col<int> &g2);

double MI(arma::Col<int> &g1, arma::Col<int> &g2);

double nlogLik(double &pAB, double &pA , double &pB, arma::Mat<double> &m);

double fmin(double ax, double bx, double tol, double &pA, double &pB, arma::Mat<double> &matrix);

double LD(arma::Col<int> &g1, arma::Col<int> &g2, const string &method="r^2");

double methods(arma::Col<int> &g1, arma::Col<int> &g2, const string &method);

arma::Mat<double> methods_mat(arma::Mat<int> &geno012, const string &method ="RMI", int thread=1);

vector<vector<string>> calc_decay(arma::Mat<int> &geno012, vector<vector<string>> &snpinfo , const string &method, int maxdist=300);

vector<vector<string>> circos(arma::Mat<double> &m, vector<vector<string>> &snpinfo , const string &method, double threshold=0.2);

void by_chrom(arma::Mat<double> &m, vector<vector<string>> &chrinfo);