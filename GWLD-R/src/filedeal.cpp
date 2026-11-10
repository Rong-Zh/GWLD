#define ARMA_64BIT_WORD 1
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../inst/include/gwld/vcf.hpp"
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// Validate ordinary genomic positions without changing the R column type.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector parse_positions(Rcpp::CharacterVector values) {
    Rcpp::NumericVector result(values.size());
    for (R_xlen_t i=0;i<values.size();++i) {
        if (Rcpp::CharacterVector::is_na(values[i])) Rcpp::stop("Site positions must not be missing");
        result[i]=static_cast<double>(GWLD::position(Rcpp::as<std::string>(values[i])));
        if (i%4096==0) Rcpp::checkUserInterrupt();
    }
    return result;
}

// [[Rcpp::export]]
arma::Mat<int> read_bed(std::string bedfile, int nSample)
{
    if (nSample <= 0) Rcpp::stop("nSample must be positive");
    std::ifstream infile(bedfile, std::ios::in | std::ios::binary);
    if (!infile)
    {
        throw Rcpp::exception("Not find *.bed file");
    }
    unsigned char magic[3];
    infile.read(reinterpret_cast<char*>(magic), 3);
    if (!infile || magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01)
        Rcpp::stop("BED must be PLINK SNP-major format");
    infile.seekg(0, std::ifstream::end);
    const std::streamoff len = infile.tellg() - std::streamoff(3);
    const std::size_t bytes = (static_cast<std::size_t>(nSample) + 3) / 4;
    if (len < 0 || len % bytes) Rcpp::stop("Truncated BED file");
    const arma::uword nSnps = len / bytes;
    infile.seekg(3, std::ifstream::beg);
    std::vector<unsigned char> buffer(bytes);
    arma::Mat<int> geno012(nSample, nSnps);
    const int decode[4] = {2, NA_INTEGER, 1, 0};
    for (arma::uword j = 0; j < nSnps; ++j) {
        infile.read(reinterpret_cast<char*>(buffer.data()), bytes);
        if (!infile) Rcpp::stop("Truncated BED file");
        int* target = geno012.colptr(j);
        for (int i = 0; i < nSample; ++i)
            target[i] = decode[(buffer[i/4] >> (2*(i%4))) & 3];
    }
    return geno012;
}

// [[Rcpp::export]]
Rcpp::List read_vcf_matrix(std::string filename,int threads=1) {
    auto data=GWLD::read_vcf_data(filename,threads);
    if (data.sites.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        Rcpp::stop("Too many variants for an R matrix");
    Rcpp::IntegerMatrix genotypes(data.samples.size(),data.sites.size());
    std::transform(data.genotypes.begin(),data.genotypes.end(),genotypes.begin(),
        [](std::int8_t code) { return code<0 ? NA_INTEGER : static_cast<int>(code); });
    Rcpp::CharacterVector chr(data.sites.size()),ids(data.sites.size());
    Rcpp::NumericVector positions(data.sites.size());
    for (std::size_t i=0;i<data.sites.size();++i) {
        chr[i]=data.sites[i].chr; ids[i]=data.sites[i].id; positions[i]=data.sites[i].pos;
    }
    Rcpp::rownames(genotypes)=Rcpp::wrap(data.samples); Rcpp::colnames(genotypes)=ids;
    return Rcpp::List::create(Rcpp::Named("Info")=Rcpp::DataFrame::create(
        Rcpp::Named("CHROM")=chr,Rcpp::Named("POS")=positions,Rcpp::Named("ID")=ids),
        Rcpp::Named("Genotype")=genotypes);
}

// The legacy character-table API is an output adapter, not another parser.
// [[Rcpp::export]]
Rcpp::CharacterMatrix read_vcf(std::string filename,std::string genotype) {
    if (genotype!="None" && genotype!="char" && genotype!="int" && genotype!="allele")
        Rcpp::stop("genotype must be None, char, int or allele");
    GWLD::VcfReader reader(filename);
    const auto samples=reader.sample_names();
    const bool paired=genotype=="allele";
    const std::size_t columns=9+samples.size()*(paired ? 2 : 1);
    if (columns>static_cast<std::size_t>(std::numeric_limits<int>::max())) Rcpp::stop("Too many output columns");
    std::vector<std::vector<std::string>> rows;
    while (reader.next(genotype!="None")) {
        auto fields=reader.text_fields();
        if (genotype!="None") {
            const auto alleles=reader.alleles();
            fields.resize(columns);
            // Write backwards because allele mode expands each sample to two columns.
            for (std::size_t i=samples.size();i-->0;) {
                const int code=reader.genotypes()[i];
                if (paired) {
                    const int a=reader.allele(i,0),b=reader.allele(i,1);
                    fields[9+2*i]=a<0 ? "NA" : std::to_string(a);
                    fields[10+2*i]=b<0 ? "NA" : std::to_string(b);
                } else if (code<0) fields[9+i]="NA";
                else if (genotype=="int") fields[9+i]=std::string(1,static_cast<char>('0'+code));
                else fields[9+i]=alleles[reader.allele(i,0)]+"/"+alleles[reader.allele(i,1)];
            }
        }
        rows.push_back(std::move(fields));
        if (rows.size()%4096==0) Rcpp::checkUserInterrupt();
    }
    if (rows.size()>static_cast<std::size_t>(std::numeric_limits<int>::max())) Rcpp::stop("Too many output rows");
    Rcpp::CharacterMatrix result(rows.size(),columns);
    const char* fixed[]={"CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"};
    Rcpp::CharacterVector names(columns);
    for (int j=0;j<9;++j) names[j]=fixed[j];
    for (std::size_t i=0;i<samples.size();++i) {
        if (paired) { names[9+2*i]=samples[i]+"_A"; names[10+2*i]=samples[i]+"_a"; }
        else names[9+i]=samples[i];
    }
    for (std::size_t i=0;i<rows.size();++i) for (std::size_t j=0;j<columns;++j) {
        if (rows[i][j]=="NA") result(i,j)=NA_STRING;
        else result(i,j)=rows[i][j];
    }
    Rcpp::colnames(result)=names;
    return result;
}
