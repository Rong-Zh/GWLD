// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>

using namespace Rcpp;
using namespace arma;

void bed_decode(char *out, const char *in, const int n) {
    char tmp, geno;
    int a1, a2;
    for(int i = 0 ; i < n ; ++i) {
        tmp = in[i];
        int k = 4 * i;
        /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
        * allele 2. The final genotype is the sum of the alleles, except for 01
        * which denotes missing. 3 for NA*/
        geno = (tmp & 3);
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        out[k] = (geno == 1) ? 3 : a1 + a2;
        k++;

        geno = (tmp & 12) >> 2; 
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        out[k] = (geno == 1) ? 3 : a1 + a2;
        k++;

        geno = (tmp & 48) >> 4; 
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        out[k] = (geno == 1) ? 3 : a1 + a2;
        k++;
        
        geno = (tmp & 192) >> 6; 
        a1 = !(geno & 1);
        a2 = !(geno >> 1);
        out[k] = (geno == 1) ? 3 : a1 + a2;
   }
}


// [[Rcpp::export]]
arma::Mat<int> read_bed(std::string bedfile, int nSample) {
    std::ifstream infile(bedfile, std::ios::in | std::ios::binary);
    if(!infile) {
        throw Rcpp::exception("Not find *.bed file");
    }
    infile.seekg(0, std::ifstream::end);
    // 从第三个字节开始
    long int len = (int)infile.tellg() - 3;
    // 向上取整，计算有多少个SNP
    int np = (int)ceil((double)nSample/4);
    int nSnps = len / np;
    infile.seekg(3, std::ifstream::beg);
    char* tmp = new char[np];
    char* tmp2 = new char[np * 4];
    arma::Mat<int> geno012(nSample, nSnps);
    arma::Col<int> tmp3(nSample);
    //用3来表示缺失,二进制为11
    for(int j = 0 ; j < nSnps ; j++) {
            // 读取二进制的基因型
            infile.read((char*)tmp, sizeof(char) * np);
            // 解码
            bed_decode(tmp2, tmp, np);
            for(int i = 0 ; i < nSample ; i++) {
                tmp3(i) = (int)tmp2[i];      
            }          
            geno012.col(j) = tmp3;
    }
    infile.close();//关闭
    geno012.replace(3, NA_REAL); //用-1 表示缺失
    return geno012;
}