#include "legacy/utils.h"
#include "legacy/utils_reference.h"
#include "dataset.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <cstdio>

void require(bool value,const char* message) { if(!value) throw std::runtime_error(message); }
std::string contents(const std::string& path) { std::ifstream in(path,std::ios::binary); return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()}; }
void same(const std::vector<double>& a,const std::vector<double>& b) {
    require(a.size()==b.size(),"Frequency shape differs");
    for(std::size_t i=0;i<a.size();++i) require((std::isnan(a[i])&&std::isnan(b[i])) || a[i]==b[i],"Frequency differs");
}
int main(int argc,char** argv) {
 try {
    const std::string prefix=argc>1 ? argv[1] : "utils-test";
    std::mt19937 random(2026);
    for(const std::string line:{"",",",",a,,b,","a\tb c\t","one","\t\t\t"}) {
        std::vector<std::string> old={"existing"},fresh=old;
        legacy_reference::split(old,line,{'\t',',',' '}); split(fresh,line,{'\t',',',' '});
        require(old==fresh,"Split tokens differ");
    }
    std::vector<std::string> alias={"a,b"}; split(alias,alias[0],{','});
    require(alias==std::vector<std::string>{"a,b","a","b"},"Aliased split source was invalidated");
    std::vector<std::string> no_delimiters; split(no_delimiters,"whole",{});
    require(no_delimiters==std::vector<std::string>{"whole"},"Empty delimiters failed");
    for(unsigned trial=0;trial<1000;++trial) {
        std::string line; for(unsigned i=0;i<100;++i) line+="ACGT \t\r\n\f\v"[random()%10];
        std::istringstream in(line); std::vector<std::string> expected; std::string token;
        while(in>>token) expected.push_back(token);
        require(expected==GWLD::fields(line),"Production whitespace parsing differs");
    }
    arma::Mat<int> genotypes(1001,13);
    for(auto& x:genotypes) x=static_cast<int>(random()%4)-1;
    auto old_freq=legacy_reference::stats_freq(genotypes),fresh_freq=stats_freq(genotypes);
    for(std::size_t i=0;i<old_freq.size();++i) same(old_freq[i],fresh_freq[i]);
    arma::Col<int> a=genotypes.col(0),b=genotypes.col(1);
    require(arma::approx_equal(table(a,b,3,3),legacy_reference::table(a,b,3,3),"absdiff",0),"Contingency table differs");
    a.fill(-1); same(allele_freq(a),{NAN,NAN,1});
    a[0]=std::numeric_limits<int>::min(); same(allele_freq(a),{NAN,NAN,1});
    a[0]=3; bool failed=false; try{allele_freq(a);}catch(const std::invalid_argument&){failed=true;} require(failed,"Invalid genotype accepted");
    require(getchrinfo({{"1"},{"1"},{"2"}})==std::vector<std::vector<std::string>>{{"1","0","1","2"},{"2","2","2","1"}},"Last chromosome was lost");
    require(getchrinfo({}).empty(),"Empty chromosome metadata failed");
    require(getsnpinfo({{"1","10","id"}},{2,0})==std::vector<std::vector<std::string>>{{"id","1"}},"Metadata selection failed");
    std::string input=prefix+"-input.tsv",old_path=prefix+"-old.tsv",new_path=prefix+"-new.tsv";
    {std::ofstream out(input);out<<"\na\t\tb\t\n\t\nlast\titem\r\n";}
    require(readfile(input,'\t')==legacy_reference::readfile(input,'\t'),"File fields differ");
    arma::Mat<double> numbers(1001,17);
    for(auto& x:numbers) x=std::ldexp(static_cast<double>(random())/random.max(),static_cast<int>(random()%2000)-1000);
    numbers[0]=NAN;numbers[1]=-NAN;numbers[2]=INFINITY;numbers[3]=-INFINITY;numbers[4]=-0.0;
    legacy_reference::writearma(numbers,old_path,'\t');writearma(numbers,new_path,'\t');
    require(contents(old_path)==contents(new_path),"Numeric file formatting differs");
    std::vector<std::vector<std::string>> text={{"1","position","id"},{"2","20","name"}};
    legacy_reference::writefile(text,old_path,'\t');writefile(text,new_path,'\t');
    require(contents(old_path)==contents(new_path),"Text table formatting differs");
    std::vector<std::vector<double>> numeric={{0,-0.0,NAN},{1.25,INFINITY,-INFINITY}};
    legacy_reference::writefile(numeric,old_path,'\t');writefile(numeric,new_path,'\t');
    require(contents(old_path)==contents(new_path),"Numeric table formatting differs");
    failed=false;try{readfile(prefix+"-nonexistent",'\t');}catch(const std::runtime_error&){failed=true;}
    require(failed,"Missing input must throw instead of terminating the process");
    std::vector<std::vector<std::string>> empty;writefile(empty,new_path,'\t'); require(contents(new_path).empty(),"Empty output failed");
    std::vector<std::string> alleles={"G","A","G","A","0"};
    require(counts_allele(alleles)==std::vector<std::pair<std::string,int>>{{"A",2},{"G",2},{"0",1}},"Allele ties are not deterministic");
    for(const auto& path:{input,old_path,new_path}) std::remove(path.c_str());
    std::cout<<"Utility parsing, frequencies, tables, file bytes and metadata boundaries passed\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
