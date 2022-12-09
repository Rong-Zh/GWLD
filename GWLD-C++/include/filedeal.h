#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <vector>
#include <string>

using namespace std;
using namespace arma;

void decode_plink(char *out, const char *in, const int n);

arma::Mat<int> plink2arma(const string &bedfile, int nSample);

vector<vector<string>> plink2vcf(vector<vector<string>> &bim_array, vector<vector<string>> &fam_array, arma::Mat<int> &geno012, bool code012=false);

arma::Mat<int> ped2arma(vector<vector<string>> &ped_array);

void pedmap2vcf(vector<vector<string>> &ped_array, vector<vector<string>> &map_array, vector<vector<string>> &vcf, vector<vector<double>> &alleles_freq);

void read_vcf(vector<vector<string>> &strArray, vector<vector<string>> &allele_mat, vector<vector<double>> &alleles_freq, const string &vcffile);

arma::Mat<int> vcf2arma(vector<vector<string>> strArray);