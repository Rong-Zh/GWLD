#pragma once
#include <armadillo>
#include <string>
#include <utility>
#include <vector>

std::vector<char> is_any_of(const std::string& str);
void split(std::vector<std::string>& result,const std::string& str,const std::vector<char>& delimiters);
std::vector<std::vector<std::string>> readfile(const std::string& filename,char sep);
template<class T> void writefile(const T& data,const std::string& filename,char sep);
void writearma(const arma::Mat<double>& data,const std::string& filename,char sep);
bool cmp(const std::pair<std::string,int>& a,const std::pair<std::string,int>& b);
std::vector<std::pair<std::string,int>> counts_allele(const std::vector<std::string>& alleles);
std::vector<double> allele_freq(const arma::Col<int>& genotype);
arma::Mat<double> table(const arma::Col<int>& a,const arma::Col<int>& b,int rows,int cols);
std::vector<std::vector<double>> stats_freq(const arma::Mat<int>& genotypes);
std::vector<std::vector<std::string>> getsnpinfo(const std::vector<std::vector<std::string>>& data,const std::vector<int>& columns);
// Return contiguous chromosome groups as name, first index, last index and count.
std::vector<std::vector<std::string>> getchrinfo(const std::vector<std::vector<std::string>>& sites);
