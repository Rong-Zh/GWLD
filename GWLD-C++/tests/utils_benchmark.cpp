#include "legacy/utils.h"
#include "legacy/utils_reference.h"
#include "dataset.hpp"
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <cstdio>

std::string bytes(const std::string& path) {std::ifstream in(path,std::ios::binary);return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()};}
std::vector<std::string> reference_fields(const std::string& line) {
    std::istringstream in(line);std::string word;std::vector<std::string> result;
    while(in>>word) result.push_back(word);
    return result;
}
void measure(const std::string& name,const std::function<void()>& before,const std::function<void()>& after,const std::function<void()>& check) {
    std::vector<double> old_times,new_times;
    for(unsigned repeat=0;repeat<4;++repeat) {
        const auto run=[&](bool old) {
            const auto start=std::chrono::steady_clock::now();
            (old ? before:after)();
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            if(repeat) (old ? old_times:new_times).push_back(seconds);
        };
        run(repeat%2==0);run(repeat%2!=0);check();
    }
    std::sort(old_times.begin(),old_times.end());std::sort(new_times.begin(),new_times.end());
    std::cout<<name<<'\t'<<old_times[1]<<'\t'<<new_times[1]<<'\t'<<old_times[1]/new_times[1]<<std::endl;
}
int main(int argc,char** argv) {
 try {
    const std::string prefix=argc>1 ? argv[1]:"utils-perf";
    std::string input=prefix+"-input.tsv",old_path=prefix+"-old.tsv",new_path=prefix+"-new.tsv";
    std::string line="1,1200042;rs1234\tA G,0\t.";
    std::vector<char> delimiters={' ', '\t', ',', ';'};
    std::vector<std::string> old_tokens,new_tokens;
    std::cout<<"operation\tbefore_seconds\tafter_seconds\tspeedup\n";
    measure("split",[&](){for(int i=0;i<200000;++i){old_tokens.clear();legacy_reference::split(old_tokens,line,delimiters);}},
        [&](){for(int i=0;i<200000;++i){new_tokens.clear();split(new_tokens,line,delimiters);}},
        [&](){if(old_tokens!=new_tokens)throw std::runtime_error("Split mismatch");});
    std::string ped_line="family sample 0 0 1 -9";
    for(int i=0;i<2000;++i) ped_line+=" A G";
    measure("production_fields",[&](){for(int i=0;i<200;++i)old_tokens=reference_fields(ped_line);},
        [&](){for(int i=0;i<200;++i)new_tokens=GWLD::fields(ped_line);},
        [&](){if(old_tokens!=new_tokens)throw std::runtime_error("Fields mismatch");});
    {std::ofstream out(input);for(int row=0;row<40000;++row)out<<"1\trs1234\t0\t10000\tA\tG\tA\tA\tG\tG\t0\t0\n";}
    std::vector<std::vector<std::string>> old_rows,new_rows;
    measure("readfile",[&](){old_rows=legacy_reference::readfile(input,'\t');},[&](){new_rows=readfile(input,'\t');},
        [&](){if(old_rows!=new_rows)throw std::runtime_error("Read mismatch");});
    arma::Mat<int> genotypes(2000,2000);std::mt19937 random(2026);
    for(auto& x:genotypes)x=static_cast<int>(random()%4)-1;
    std::vector<std::vector<double>> old_freq,new_freq;
    measure("stats_freq",[&](){old_freq=legacy_reference::stats_freq(genotypes);},[&](){new_freq=stats_freq(genotypes);},
        [&](){if(old_freq!=new_freq)throw std::runtime_error("Frequency mismatch");});
    std::vector<std::string> alleles(500000);
    for(auto& x:alleles){const auto code=random()%100;x=code<5 ? "0":code<70 ? "A":"G";}
    std::vector<std::pair<std::string,int>> old_counts,new_counts;
    measure("counts_allele",[&](){old_counts=legacy_reference::counts_allele(alleles);},[&](){new_counts=counts_allele(alleles);},
        [&](){if(old_counts!=new_counts)throw std::runtime_error("Allele count mismatch");});
    arma::Mat<double> numbers(1024,1024);
    for(auto& x:numbers)x=random()%2 ? static_cast<double>(random())/random.max():NAN;
    measure("writearma",[&](){legacy_reference::writearma(numbers,old_path,'\t');},[&](){writearma(numbers,new_path,'\t');},
        [&](){if(bytes(old_path)!=bytes(new_path))throw std::runtime_error("Write mismatch");});
    for(const auto& path:{input,old_path,new_path})std::remove(path.c_str());
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
