#pragma once
#include <sstream>
#include "dataset.hpp"
namespace GWLD { namespace legacy {
inline void finalize(Dataset& data,const std::vector<int>& flat) {
    if (data.sites.empty() || data.samples.empty()) throw std::runtime_error("Input has no variants or samples");
    data.genotypes.set_size(data.samples.size(),data.sites.size());
    std::copy(flat.begin(),flat.end(),data.genotypes.memptr());
}
inline std::vector<std::string> fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> result;
    std::string word;
    while (stream>>word) result.push_back(word);
    return result;
}
inline Dataset load_vcf(const std::string& path) {
    auto stream=input(path);
    Dataset data;
    std::vector<int> flat;
    std::string line;
    std::size_t number=0;
    while (std::getline(stream,line)) {
        ++number;
        if (line.empty() || line.compare(0,2,"##")==0) continue;
        auto row=fields(line);
        if (line.compare(0,6,"#CHROM")==0) {
            if (!data.samples.empty() || row.size()<10) throw std::runtime_error("Invalid VCF sample header");
            data.samples.assign(row.begin()+9,row.end()); continue;
        }
        if (line[0]=='#') continue;
        const auto error="Invalid VCF record at line "+std::to_string(number);
        if (data.samples.empty() || row.size()!=data.samples.size()+9) throw std::runtime_error(error);
        if (row[4].find(',')!=std::string::npos) throw std::runtime_error(error+": only biallelic variants are supported");
        std::vector<std::string> format;
        std::istringstream formats(row[8]); std::string field;
        while (std::getline(formats,field,':')) format.push_back(field);
        const auto gt=std::find(format.begin(),format.end(),"GT");
        if (gt==format.end()) throw std::runtime_error(error+": FORMAT has no GT");
        const auto gtindex=static_cast<std::size_t>(gt-format.begin());
        data.sites.push_back({row[0],row[2],position(row[1])});
        data.alleles.push_back({{row[3],row[4]}});
        for (std::size_t sample=9; sample<row.size(); ++sample) {
            std::istringstream values(row[sample]); std::string genotype;
            for (std::size_t k=0; k<=gtindex; ++k)
                if (!std::getline(values,genotype,':')) throw std::runtime_error(error+": incomplete FORMAT");
            int code=-1;
            if (genotype=="0/0" || genotype=="0|0") code=0;
            else if (genotype=="0/1" || genotype=="1/0" || genotype=="0|1" || genotype=="1|0") code=1;
            else if (genotype=="1/1" || genotype=="1|1") code=2;
            else if (genotype!="." && genotype!="./." && genotype!=".|."
                     && genotype!="0/." && genotype!="./0" && genotype!="1/." && genotype!="./1"
                     && genotype!="0|." && genotype!=".|0" && genotype!="1|." && genotype!=".|1")
                throw std::runtime_error(error+": unsupported diploid GT "+genotype);
            flat.push_back(code);
        }
    }
    if (stream.bad()) throw std::runtime_error("Error reading VCF: "+path);
    finalize(data,flat); return data;
}
} } // namespace GWLD::legacy
