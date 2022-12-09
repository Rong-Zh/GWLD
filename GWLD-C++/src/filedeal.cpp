#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <armadillo>
#include <float.h>
#include "utils.h"

using namespace std;
using namespace arma;

void decode_plink(char *out, const char *in, const int n) {
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

arma::Mat<int> plink2arma(const string &bedfile, int nSample) {
    std::ifstream infile(bedfile, std::ios::in | std::ios::binary);
    if(!infile) {
        std::cerr << bedfile <<".bed 读取失败" << std::endl;
        throw std::runtime_error("IO Error");
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
            decode_plink(tmp2, tmp, np);
            for(int i = 0 ; i < nSample ; i++) {
                tmp3(i) = (int)tmp2[i];      
            }          
            geno012.col(j) = tmp3;
    }
    infile.close();//关闭
    geno012.replace(3, -1); //用-1 表示缺失
    return geno012;
}

vector<vector<string>> plink2vcf(vector<vector<string>> &bim_array, vector<vector<string>> &fam_array, arma::Mat<int> &geno012, bool code012=false) {
    int nSNP = bim_array.size();
    int nSample = fam_array.size();
    vector<vector<string>> vcf;
    vector<string> vcf_vector(nSample+9);
    for(int i=0; i<nSNP; i++) {
        vcf_vector[0] = bim_array[i][0];
        vcf_vector[1] = bim_array[i][3];
        vcf_vector[2] = bim_array[i][1];
        vcf_vector[3] = bim_array[i][4];
        vcf_vector[4] = bim_array[i][5];
        vcf_vector[5] = ".";
        vcf_vector[6] = ".";
        vcf_vector[7] = "PR";
        vcf_vector[8] = "GT";
        for(int j=0; j<nSample; j++) { 
            if (!code012 && geno012(j,i)== 0) {
                vcf_vector[j+9] = "0/0";
            } else if(!code012 && geno012(j,i)== 1) {
                vcf_vector[j+9] = "0/1";
            } else if(!code012 && geno012(j,i)== 2) {
                vcf_vector[j+9] = "1/1";
            } else if(!code012 && geno012(j,i)== -1) {
                vcf_vector[j+9] = "./.";
            } else {
                vcf_vector[j+9] = to_string(geno012(j,i));
            }
        }
        vcf.push_back(vcf_vector);
    }
    return vcf;
}

arma::Mat<int> ped2arma(vector<vector<string>> &ped_array) {
    //查看字符数组的维度
    int nSample = ped_array.size();
    long int nSNP = (ped_array.front().size()-6)/2;
    //前6列,FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype. 总共nSNP*2+6列
    vector<string> alleles(2*nSample); //存放单个等位基因
    vector<string> genotype(nSample);
    arma::Mat<int> geno012(nSample, nSNP);
    for(int i=0; i<nSNP; i++) {
        for(int j=0; j<nSample; j++){
            alleles[2*j] = ped_array[j][2*i+6]; //第一个等位基因
            alleles[2*j+1] = ped_array[j][2*i+7]; //第二个等位基因
            genotype[j] = ped_array[j][2*i+6] + ped_array[j][2*i+7];//把两个等位基因和成一个基因型
        }
        vector<pair<string,int>>  allele_sorted = counts_allele(alleles);
        string A1 = allele_sorted[0].first;//major allele
        string A2 = allele_sorted[1].first;//major allele
        //存入数组中
        //用0, 1, 2, -1 替换碱基基因型
        for(int index =0; index<genotype.size(); index++) { //行填充
            if(genotype[index]==(A1+A1)) {
                geno012(index, i) = 0;
            } else if(genotype[index]==(A1+A2) || genotype[index]==(A2+A1)) {
                geno012(index, i) = 1;
            } else if(genotype[index]==(A2+A2)) {
                geno012(index, i) = 2;
            } else {
                geno012(index, i) = -1;//-1 表示缺失
            }
        }
    }
    return geno012;
}

void pedmap2vcf(vector<vector<string>> &ped_array, vector<vector<string>> &map_array, vector<vector<string>> &vcf, vector<vector<double>> &alleles_freq) { 
    // ped_array    输入文件ped
    // map_array    输入文件map
    // vcf          输出 是转换基因型后的数组
    // alleles_freq 输出 存放等位基因频率结果
    //查看字符数组的维度
    int nSample = ped_array.size();
    long int nSNP = (ped_array.front().size()-6)/2;
    long int map_nSNP = map_array.size(); // 查看行数
    //前6列,FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype. 总共2*nSNP+6列
    vector<string> alleles(2*nSample);    // 存放单个等位基因
    vector<double> alleles_vector(3);     // 存放等位基因major allele, minor allele, NA的频率
    vector<string> genotype(nSample);     // 存放每个位点的基因型, 有nSample个         // 存放转换后的基因型
    vector<string> vcf_vector(nSample+9);
    string A1, A2;
    // 构造vcf文件中的每一行
    //判断ped和map文件是否对应
    if(nSNP == map_nSNP) {
        for(int i=0; i<nSNP; i++) { 
            for(int j=0; j<nSample; j++){
                alleles[2*j] = ped_array[j][2*i+6]; //第一个等位基因
                alleles[2*j+1] = ped_array[j][2*i+7]; //第二个等位基因
                genotype[j] = ped_array[j][2*i+6] + ped_array[j][2*i+7];//把两个等位基因和成一个基因型
            }
            vector<pair<string,int>>  tmp = counts_allele(alleles);
            A1 = tmp[0].first;//major allele
            A2 = tmp[1].first;//minor allele
            int sum = tmp[0].second+tmp[1].second;//计算两种等位基因总和
            alleles_vector[0] = (double)tmp[0].second/(sum);//计算major allele 频率
            alleles_vector[1] = (double)tmp[1].second/(sum);//计算minor allele 频率
            tmp.size()==2 ? alleles_vector[2] = 0.0 : alleles_vector[2] = (double)tmp[2].second/(2*nSample);
            //存入数组中
            alleles_freq.push_back(alleles_vector);
            vcf_vector[0] = map_array[i][0];
            vcf_vector[1] = map_array[i][3];
            vcf_vector[2] = map_array[i][1];
            vcf_vector[3] = A1;
            vcf_vector[4] = A2;
            vcf_vector[5] = ".";
            vcf_vector[6] = ".";
            vcf_vector[7] = "PR";
            vcf_vector[8] = "GT";
            //转换成vcf文件
            for(int index=0; index<genotype.size(); index++) { //行填充
                if(genotype[index]==(A1+A1)) {
                    vcf_vector[index+9] = "0/0";
                } else if(genotype[index]==(A1+A2) || genotype[index]==(A2+A1)) {
                    vcf_vector[index+9] = "0/1";
                } else if(genotype[index]==(A2+A2)) {
                    vcf_vector[index+9] = "1/1";
                } else {
                    vcf_vector[index+9] = "./."; //vcf中用./.表示缺失
                }
            }
            vcf.push_back(vcf_vector);
        }
    } else {
        std::cerr << "ped文件中的SNP和map文件中SNP数不相等" << endl;
        return;
    }
}

void read_vcf(vector<vector<string>> &strArray, vector<vector<string>> &allele_mat, vector<vector<double>> &alleles_freq, const string &vcffile) {
    //读取文件
    ifstream infile;
    infile.open(vcffile, ios::in);
    if (!infile.is_open())
	{
		std::cerr << "读取文件失败" << endl;
        return;
	}
    string lineStr;
    string A1, A2; 
    vector<double> alleles_vector(3); //转换为-1,0,1,2基因型,-1为缺失基因型
	while (getline(infile, lineStr))
	{
        if(lineStr[0]=='#') //跳过注释的行
            continue;
        stringstream ss(lineStr);
        string str;
        vector<string> lineArray;
        while (getline(ss, str, '\t')) {
            lineArray.push_back(str);
        }
        vector<string> line;
        split(line, lineStr, is_any_of("\t/|"));
        vector<string> gt(line.begin()+9, line.end());
        int nSample = line.size()-9;
        replace(gt.begin(), gt.end(), string("."), string("-1"));//把.替换成-1
        vector<pair<string,int>>  tmp = counts_allele(gt);
        A1 = tmp[0].first;//major allele
        A2 = tmp[1].first;//minor allele
        int sum = tmp[0].second + tmp[1].second;//等位基因的总个数
        alleles_vector[0] = (double)tmp[0].second/(sum);//计算major allele 频率
        alleles_vector[1] = (double)tmp[1].second/(sum);//计算minor allele 频率
        tmp.size()==2 ? alleles_vector[2] = 0.0 : alleles_vector[2] = (double)tmp[2].second/(2*nSample);

        for(int j=9; j<lineArray.size(); j++) {
            if(lineArray[j]==(A1+"/"+A1)|| lineArray[j]== (A1+"|"+A1)) {
                lineArray[j] = "0";
            } else if (lineArray[j]==(A1+"/"+A2) || lineArray[j]==(A2+"/"+A1) || lineArray[j]==(A1+"|"+A2) || lineArray[j]==(A2+"|"+A1)) {
                lineArray[j] = "1";
            } else if (lineArray[j]==(A2+"/"+A2) || lineArray[j]==(A2+"|"+A2)) {
                lineArray[j] = "2";
            } else {
                lineArray[j] = "-1";
            }
        }
        alleles_freq.push_back(alleles_vector);
        allele_mat.push_back(gt);
        strArray.push_back(lineArray);
    }
}

arma::Mat<int> vcf2arma(vector<vector<string>> strArray) {
    int nSNP = strArray.size();
    int rows = strArray.front().size();
    arma::Mat<int> geno012(rows-9, nSNP);//行是位点, 列为样本
    for(int i=9; i<rows; i++) {
        for(int j=0; j<nSNP; j++) {
                geno012(i-9, j) = std::stoi(strArray[j][i]);
        }
    }
    return geno012;
}