// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

using namespace Rcpp;
using namespace arma;
using namespace std;


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
    // 用3来表示缺失,二进制为11
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

std::vector<std::vector<std::string>> read_table(std::string filename, std::string genotype) {

    ifstream infile;
    infile.open(filename, ios::in);
    if (!infile.is_open()) {
        throw Rcpp::exception("*.vcf file loading fail");
	}
    std::string lineStr;
    std::vector<std::vector<std::string>> strArray;
    while (getline(infile, lineStr))
	{   // 跳过注释的行
        if(lineStr.find("##")!=std::string::npos) 
            continue;
        stringstream ss(lineStr);
        std::string str;
        std::vector<std::string> lineArray;
        std::vector<std::string> row;
        while (getline(ss, str, '\t')) {
            lineArray.push_back(str);
        }
        if(genotype=="None") {
            row = lineArray;
        } else if(genotype=="char") {
            for(size_t i=0; i<lineArray.size(); i++) {
                if(lineArray[i]=="0/0" || lineArray[i]=="0|0") {
                    row.push_back(lineArray[3]+"/"+lineArray[3]);
                } else if(lineArray[i]=="0/1" || lineArray[i]=="0|1") {
                    row.push_back(lineArray[3]+"/"+lineArray[4]);
                } else if(lineArray[i]=="1/1" || lineArray[i]=="1|1") {
                    row.push_back(lineArray[4]+"/"+lineArray[4]);
                } else if(lineArray[i]=="./." || lineArray[i]==".|.") {
                    row.push_back("NA");
                } else {
                    row.push_back(lineArray[i]);
                }

            }
        } else if(genotype=="int") {
            for(size_t i=0; i<lineArray.size(); i++) {
                if(lineArray[i]=="0/0" || lineArray[i]=="0|0") {
                    row.push_back("0");
                } else if(lineArray[i]=="0/1" || lineArray[i]=="0|1") {
                    row.push_back("1");
                } else if(lineArray[i]=="1/1" || lineArray[i]=="1|1") {
                    row.push_back("2");
                } else if(lineArray[i]=="./." || lineArray[i]==".|.") {
                    row.push_back("NA");
                } else {
                    row.push_back(lineArray[i]);
                }  
            }
        } else if(genotype=="allele") {
            for(size_t i=0; i<lineArray.size(); i++) {
                if(lineArray[i]=="0/0" || lineArray[i]=="0|0") {
                    row.push_back("0");
                    row.push_back("0");
                } else if(lineArray[i]=="0/1" || lineArray[i]=="0|1") {
                    row.push_back("0");
                    row.push_back("1");
                } else if(lineArray[i]=="1/1" || lineArray[i]=="1|1") {
                    row.push_back("1");
                    row.push_back("1");
                } else if(lineArray[i]=="./." || lineArray[i]==".|.") {
                    row.push_back("NA");
                    row.push_back("NA");
                } else if(i>=9) {
                    row.push_back(lineArray[i]+"_A");
                    row.push_back(lineArray[i]+"_a");
                } else {
                    row.push_back(lineArray[i]);
                }
            }
        } else {
            throw Rcpp::exception("Unknown genotype! Input one of methods(None, char, int, allele) ");
        }
        strArray.push_back(row);
    }
    infile.close();
    return strArray;
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix read_vcf(std::string filename, std::string genotype) {
    std::vector<std::vector<std::string>> cMat = read_table(filename, genotype);
    Rcpp::CharacterVector ColName(cMat[0].begin(), cMat[0].end());
    if (cMat[0][0].find('#')!=std::string::npos) {
        //去掉首行行名中的#注释
        ColName[0] = cMat[0][0].replace(cMat[0][0].find('#'), 1, "");
    }
    int nrows = cMat.size(), ncols=cMat.front().size();
    Rcpp::CharacterMatrix rMat(nrows-1, ncols);
    //把c++矩阵转换为R语言矩阵
    for(int row=1; row<nrows; row++ ){
        for(int col=0; col<ncols; col++) {
            if (cMat[row][col]=="NA") {
                //把字符NA替换成R语言NA缺失值
                rMat(row-1, col) = NA_STRING; 
            } else {
                rMat(row-1, col) = cMat[row][col];
            }
        }
    }
    //设置列名  
    colnames(rMat) = ColName;
    return rMat;
}