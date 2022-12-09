#include <iostream>
#include <string>
#include <armadillo>
#include <vector>
#include "cmdline.h"
#include "filedeal.h"
#include "utils.h"
#include "methods.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]) {
    cmdline::parser a;
    a.add<string>("bfile", 'b', "Input plink format (*.bed *.bim and *.fam) file's prefix", false);
    a.add<string>("vcf", 'v', "Input vcf format file name", false);
    a.add<string>("file", 'f', "Input plink format (*.ped and *.map) file's prefix", false);
    a.add<string>("method", 'm', "Calculation method between two sites", false, "RMI", cmdline::oneof<string>("RMI", "MI", "r2", "D", "D'"));
    a.add("by-chrom", '\0', "Calculated by chromosome");
    a.add("decay", 'd', "Output decay according to the chosen method");
    a.add<int>("max-dist", '\0', "Calculates the maximum distance between two sites. default 300kb", false, 300);
    a.add("circos", 'c', "Calculate between different chromosomes");
    a.add("code012", '\0', "Recode genotype with -1, 0 ,1, 2. -1 for missing value");
    a.add("allele-freq", '\0', "Calculate allele frequency per site");
    a.add("2vcf", '\0', "Convert plink format file to vcf format file");
    a.add<int>("threads", 'p', "threads number", false, 1);
    a.add<string>("out-prefix", 'o', "Output file name's prefix", false, "gwld");  
    a.add("help", 'h', "print help message");
    a.footer("filename ...");
    a.set_program_name("GWLD");

    bool ok=a.parse(argc, argv);

    if (argc==1 || a.exist("help")){
        cerr<<a.usage();
        return 0;
    }
    
    if (!ok){
        cerr<<a.error()<<endl<<a.usage();
        return 0;
    }

    string outprefix = a.get<string>("outfileprefix"); //输出文件名前缀
    arma::Mat<int> geno012;
    vector<vector<string>> snpinfo;
    if(a.exist("bfile")) {//读取plink .bed, .fim, .bim文件
        string prefix = a.get<string>("bfile");
        string bedfile = prefix + ".bed";//二进制文件
        string bimfile = prefix + ".bim";//行数为SNP数目
        string famfile = prefix + ".bed";//行数为样本数目
        vector<vector<string>> bim_array = readfile(bimfile, '\t');
        vector<vector<string>> fam_array = readfile(famfile, ' ');
        int nSample = fam_array.size(); 
        geno012 = plink2arma(bedfile, nSample);//0,1,2,-1基因型矩阵，行为样品，列为SNP
        vector<int> index = {0,3,1};
        snpinfo = getsnpinfo(bim_array, index);//提取位点CHROM, POS, ID 

        if(a.exist("2vcf") || a.exist("allelefreq")) {
            if(a.exist("code012")) { 
                vector<vector<string>> vcf = plink2vcf(bim_array, fam_array, geno012, true);
                string outvcf = outprefix + ".vcf";
                writefile(vcf, outvcf, '\t');
            } else if(a.exist("allelefreq")) {
                vector<vector<double>> alleles_freq = stats_freq(geno012);//等位基因频率
                string outfreq = outprefix + ".freq";
                writefile(alleles_freq, outfreq, '\t');
            } else { 
                vector<vector<string>> vcf = plink2vcf(bim_array, fam_array, geno012, false);
                string outvcf = outprefix + ".vcf";
                writefile(vcf, outvcf, '\t');
            }
        }
    } else if(a.exist("file")) {//读取plink, ped和map文件
        string prefix = a.get<string>("file");
        string pedfile = prefix + ".ped";
        string mapfile = prefix + ".map";
        vector<vector<string>> ped_array = readfile(pedfile, ' ');
        vector<vector<string>> map_array = readfile(mapfile, '\t');
        geno012 = ped2arma(ped_array);
        vector<int> index = {0,3,1};
        snpinfo = getsnpinfo(map_array, index);//提取位点CHROM, POS, ID 
        if(a.exist("2vcf") || a.exist("allelefreq")) {
            vector<vector<string>> vcf;
            vector<vector<double>> alleles_freq;
            pedmap2vcf(ped_array, map_array, vcf, alleles_freq);
        if(a.exist("code012")) {
            string outvcf = outprefix + ".vcf";
            writefile(vcf, outvcf, '\t');
        } else if(a.exist("allelefreq")) {
            string outfreq = outprefix + ".freq";
            writefile(alleles_freq, outfreq, '\t');
            } 
        }
    } else if(a.exist("vcf")) {//读取vcf文件
        string vcffile = a.get<string>("vcf");
        vector<vector<string>> strArray; //基因型0,1,2的vcf文件
        vector<vector<string>> allele_mat;
        vector<vector<double>> alleles_freq;
        read_vcf(strArray, allele_mat, alleles_freq, vcffile);
        geno012= vcf2arma(strArray);
        vector<int> index = {0,1,2};
        snpinfo = getsnpinfo(strArray, index);//提取位点CHROM, POS, ID 
        if(a.exist("code012")) {
            string outvcf = outprefix+"vcf";
            writefile(strArray, outvcf, '\t');
        } else if(a.exist("allelefreq")) {
            string outfreq = outprefix + ".freq";
            writefile(alleles_freq, outfreq, '\t');
        }
    } else {
        cerr << "input file" << endl;
        return 0;
    }

    if(a.exist("method")) {
        string method = a.get<string>("method");
        int thread = a.get<int>("thread");//线程数
        arma::Mat<double> result = methods_mat(geno012, method, thread);
        string outres = outprefix + ".heatmap";
        writearma(result, outres, '\t'); //写入文件
        if(a.exist("decay")) {
            int maxdist = a.get<int>("maxdist");
            vector<vector<string>> decay_res = calc_decay(geno012, snpinfo, method, maxdist);
            string outdecay = outprefix + ".decay";
            writefile(decay_res, outdecay, '\t');
        } 
        if(a.exist("by-chrom")) {
            vector<vector<string>> chrinfo = getchrinfo(snpinfo);
            by_chrom(result, chrinfo);
        }
        if(a.exist("circos")) {
            vector<vector<string>> circos_res = circos(result, snpinfo, method, 0.1);
            string outcircos = outprefix+".circos";
            writefile(circos_res, outcircos,'\t');
        }
    }
    return 0;
}