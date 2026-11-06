#pragma once
#include <armadillo>
#include <fstream>
#include <sstream>
#include <array>
#include "gwld/vcf.hpp"

namespace GWLD {
inline std::ifstream input(const std::string& path, bool binary=false) {
    std::ifstream stream(path,binary ? std::ios::in|std::ios::binary : std::ios::in);
    if (!stream) throw std::runtime_error("Cannot open input: "+path);
    return stream;
}
inline std::vector<std::string> fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> result;
    std::string word;
    while (stream>>word) result.push_back(word);
    return result;
}
inline Dataset load_bed(const std::string& prefix) {
    Dataset data;
    auto fam=input(prefix+".fam"), bim=input(prefix+".bim");
    std::string line;
    while (std::getline(fam,line)) {
        auto row=fields(line); if (row.empty()) continue;
        if (row.size()!=6) throw std::runtime_error("FAM requires six fields");
        data.samples.push_back(row[1]);
    }
    while (std::getline(bim,line)) {
        auto row=fields(line); if (row.empty()) continue;
        if (row.size()!=6) throw std::runtime_error("BIM requires six fields");
        data.sites.push_back({row[0],row[1],position(row[3])});
        // v1 decodes allele-1 dosage: allele 2 is REF (dosage zero).
        data.alleles.push_back({{row[5],row[4]}});
    }
    if (data.sites.empty() || data.samples.empty()) throw std::runtime_error("Empty BIM/FAM");
    auto bed=input(prefix+".bed",true);
    unsigned char magic[3]; bed.read(reinterpret_cast<char*>(magic),3);
    if (!bed || magic[0]!=0x6c || magic[1]!=0x1b || magic[2]!=0x01)
        throw std::runtime_error("BED must be PLINK SNP-major format");
    data.genotypes.set_size(data.samples.size(),data.sites.size());
    const std::size_t bytes=(data.samples.size()+3)/4;
    std::vector<unsigned char> buffer(bytes);
    const int decode[4]={2,-1,1,0};
    for (std::size_t j=0; j<data.sites.size(); ++j) {
        bed.read(reinterpret_cast<char*>(buffer.data()),bytes);
        if (!bed) throw std::runtime_error("Truncated BED file");
        int* target=data.genotypes.colptr(j);
        for (std::size_t i=0; i<data.samples.size(); ++i) target[i]=decode[(buffer[i/4]>>(2*(i%4)))&3];
    }
    if (bed.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("BED/BIM variant count mismatch");
    return data;
}
inline Dataset load_ped(const std::string& prefix) {
    Dataset data;
    auto map=input(prefix+".map"), ped=input(prefix+".ped");
    std::string line;
    while (std::getline(map,line)) {
        auto row=fields(line); if (row.empty()) continue;
        if (row.size()!=4) throw std::runtime_error("MAP requires four fields");
        data.sites.push_back({row[0],row[1],position(row[3])});
    }
    std::vector<std::vector<std::string>> rows;
    while (std::getline(ped,line)) {
        auto row=fields(line); if (row.empty()) continue;
        if (row.size()!=6+2*data.sites.size()) throw std::runtime_error("PED/MAP dimensions do not match");
        data.samples.push_back(row[1]); rows.push_back(std::move(row));
    }
    if (data.sites.empty() || rows.empty()) throw std::runtime_error("Empty PED/MAP");
    data.genotypes.set_size(rows.size(),data.sites.size());
    for (std::size_t j=0; j<data.sites.size(); ++j) {
        std::map<std::string,std::size_t> counts;
        for (const auto& row: rows) for (int k=0; k<2; ++k)
            if (row[6+2*j+k]!="0") ++counts[row[6+2*j+k]];
        if (counts.size()>2) throw std::runtime_error("PED contains a multiallelic variant");
        std::vector<std::pair<std::string,std::size_t>> ordered(counts.begin(),counts.end());
        std::stable_sort(ordered.begin(),ordered.end(),[](const std::pair<std::string,std::size_t>& a,const std::pair<std::string,std::size_t>& b){return a.second>b.second;});
        const std::string ref=ordered.empty() ? "N" : ordered[0].first;
        const std::string alt=ordered.size()>1 ? ordered[1].first : ".";
        data.alleles.push_back({{ref,alt}});
        for (std::size_t i=0; i<rows.size(); ++i) {
            const auto& a=rows[i][6+2*j]; const auto& b=rows[i][7+2*j];
            data.genotypes(i,j)=(a=="0" || b=="0") ? -1 : static_cast<int>((a!=ref)+(b!=ref));
        }
    }
    return data;
}
inline std::ofstream output(const std::string& path) {
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("Cannot create output: "+path);
    stream.exceptions(std::ios::badbit|std::ios::failbit);
    stream.precision(12);
    return stream;
}
inline void write_vcf(const Dataset& data,const std::string& path) {
    auto out=output(path);
    out<<"##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    out<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sample:data.samples) out<<'\t'<<sample;
    out<<'\n';
    const char* gt[4]={"./.","0/0","0/1","1/1"};
    for (std::size_t j=0; j<data.sites.size(); ++j) {
        const auto& s=data.sites[j];
        out<<s.chr<<'\t'<<s.pos<<'\t'<<s.id<<'\t'<<data.alleles[j][0]<<'\t'<<data.alleles[j][1]<<"\t.\t.\t.\tGT";
        for (std::size_t i=0; i<data.samples.size(); ++i) {
            const int code=data.genotypes(i,j); out<<'\t';
            out<<gt[code+1];
        }
        out<<'\n';
    }
    out.close();
}
inline void write_dosage(const Dataset& data,const std::string& path) {
    auto out=output(path);
    out<<"CHROM\tPOS\tID\tREF\tALT";
    for (const auto& sample:data.samples) out<<'\t'<<sample;
    out<<'\n';
    for (std::size_t j=0; j<data.sites.size(); ++j) {
        const auto& site=data.sites[j];
        out<<site.chr<<'\t'<<site.pos<<'\t'<<site.id<<'\t'<<data.alleles[j][0]<<'\t'<<data.alleles[j][1];
        for (std::size_t i=0; i<data.samples.size(); ++i) {
            out<<'\t';
            if (data.genotypes(i,j)<0) out<<"NA";
            else out<<data.genotypes(i,j);
        }
        out<<'\n';
    }
    out.close();
}
} // namespace GWLD
