#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

using namespace std;
using namespace arma;

vector<char> is_any_of(string str);

void split(vector<string>& result, string str, vector<char> delimiters);

vector<vector<string>> readfile(const string &filename, char sep);

template <typename T>
void writefile(T &data, string &outfilename, char sep);

void writearma(arma::Mat<double> &data, string &outfilename, char sep);

bool cmp(pair<string,int> a, pair<string, int> b);

vector<pair<string,int>> counts_allele(vector<string> &allele);

vector<double> allele_freq(arma::Col<int> &g1);

arma::Mat<double> table(arma::Col<int> &g1, arma::Col<int> &g2, int row, int col);

vector<vector<double>> stats_freq(arma::Mat<int> &geno012);

vector<vector<string>> getsnpinfo(vector<vector<string>> &data, vector<int> &col);

vector<vector<string>> getchrinfo(vector<vector<string>> &snpinfo);