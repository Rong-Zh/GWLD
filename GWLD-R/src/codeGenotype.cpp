#define ARMA_64BIT_WORD 1
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <vector>

namespace GWLD {
// R interns character strings. Count each CHARSXP once and parse each distinct
// genotype once per column; no std::string copies are made for input cells.
class GenotypeEncoder {
    struct Allele { const char* data; std::size_t size; std::uint64_t count; };
    struct Pattern {
        int first=-1,second=-1,code=NA_INTEGER;
        std::uint64_t count=0;
    };
    bool delimiters[256]={};
    bool separated;
    std::unordered_map<SEXP,int> lookup;
    std::vector<Pattern> patterns;
    std::vector<Allele> alleles;
    R_xlen_t missing=0;

    static bool missing_token(const char* text,std::size_t size) {
        return !size || (size==1 && text[0]=='.') || (size==2 && text[0]=='N' && text[1]=='A');
    }
    static bool less(const Allele& a,const Allele& b) {
        const int order=std::memcmp(a.data,b.data,std::min(a.size,b.size));
        return order<0 || (order==0 && a.size<b.size);
    }
    int allele(const char* text,std::size_t size) {
        if (missing_token(text,size)) return -1;
        for (std::size_t i=0;i<alleles.size();++i)
            if (size==alleles[i].size && std::memcmp(text,alleles[i].data,size)==0) return static_cast<int>(i);
        if (alleles.size()>=static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("Too many distinct alleles");
        alleles.push_back({text,size,0});
        return static_cast<int>(alleles.size()-1);
    }
    Pattern parse(SEXP value) {
        const char* text=CHAR(value);
        const auto size=static_cast<std::size_t>(LENGTH(value));
        Pattern pattern;
        if (missing_token(text,size)) return pattern;
        if (!separated && size==2) {
            pattern.first=allele(text,1); pattern.second=allele(text+1,1);
        } else if (separated && size==3 && delimiters[static_cast<unsigned char>(text[1])]) {
            if (delimiters[static_cast<unsigned char>(text[0])] || delimiters[static_cast<unsigned char>(text[2])])
                throw std::invalid_argument("Expected exactly two alleles; check genotype separators");
            pattern.first=allele(text,1); pattern.second=allele(text+2,1);
        } else if (separated) {
            std::size_t split=0;
            while (split<size && !delimiters[static_cast<unsigned char>(text[split])]) ++split;
            bool valid=split>0 && split+1<size;
            for (std::size_t i=split+1;i<size && valid;++i)
                if (delimiters[static_cast<unsigned char>(text[i])]) valid=false;
            if (!valid) throw std::invalid_argument("Expected exactly two alleles; check genotype separators");
            pattern.first=allele(text,split); pattern.second=allele(text+split+1,size-split-1);
        } else {
            throw std::invalid_argument("Without separators, each diploid genotype must contain two single-byte alleles");
        }
        return pattern;
    }
public:
    explicit GenotypeEncoder(const std::string& sep) : separated(!sep.empty()) {
        for (unsigned char byte:sep) delimiters[byte]=true;
        lookup.reserve(16); patterns.reserve(16); alleles.reserve(4);
    }
    void analyze(const SEXP* source,R_xlen_t size,int* output=nullptr) {
        lookup.clear(); patterns.clear(); alleles.clear(); missing=size;
        for (R_xlen_t i=0;i<size;++i) {
            const SEXP value=source[i];
            if (value==NA_STRING) { if (output) output[i]=-1; continue; }
            auto found=lookup.find(value);
            int index;
            if (found==lookup.end()) {
                if (patterns.size()>=static_cast<std::size_t>(std::numeric_limits<int>::max()))
                    throw std::length_error("Too many distinct genotypes");
                index=static_cast<int>(patterns.size());
                patterns.push_back(parse(value)); lookup.emplace(value,index);
            } else index=found->second;
            ++patterns[index].count;
            if (output) output[i]=index;
        }
        for (const auto& pattern:patterns) {
            if (pattern.first>=0) alleles[pattern.first].count+=pattern.count;
            if (pattern.second>=0) alleles[pattern.second].count+=pattern.count;
            if (pattern.first>=0 && pattern.second>=0) missing-=static_cast<R_xlen_t>(pattern.count);
        }
        int major=-1;
        for (std::size_t i=0;i<alleles.size();++i)
            if (major<0 || alleles[i].count>alleles[major].count ||
                (alleles[i].count==alleles[major].count && less(alleles[i],alleles[major]))) major=static_cast<int>(i);
        for (auto& pattern:patterns)
            if (pattern.first>=0 && pattern.second>=0)
                pattern.code=2-(pattern.first==major)-(pattern.second==major);
        if (output) for (R_xlen_t i=0;i<size;++i)
            output[i]=output[i]<0 ? NA_INTEGER : patterns[output[i]].code;
    }
    Rcpp::List frequencies(R_xlen_t size) const {
        std::vector<std::size_t> order(alleles.size());
        for (std::size_t i=0;i<order.size();++i) order[i]=i;
        std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b) {
            return alleles[a].count!=alleles[b].count ? alleles[a].count>alleles[b].count : less(alleles[a],alleles[b]);
        });
        double total=0;
        for (const auto& allele:alleles) total+=allele.count;
        Rcpp::CharacterVector names(order.size());
        Rcpp::NumericVector counts(order.size()),frequencies(order.size());
        for (std::size_t i=0;i<order.size();++i) {
            const auto& allele=alleles[order[i]];
            names[i]=std::string(allele.data,allele.size);
            counts[i]=static_cast<double>(allele.count); frequencies[i]=counts[i]/total;
        }
        Rcpp::RObject missing_count=missing<=std::numeric_limits<int>::max()
            ? Rcpp::wrap(static_cast<int>(missing)) : Rcpp::wrap(static_cast<double>(missing));
        return Rcpp::List::create(Rcpp::Named("allele")=names,Rcpp::Named("allele counts")=counts,
            Rcpp::Named("allele freq")=frequencies,Rcpp::Named("missing counts")=missing_count,
            Rcpp::Named("missing rate")=size ? static_cast<double>(missing)/size : NA_REAL);
    }
};
} // namespace GWLD

//' Count allele frequencies in diploid character genotypes
//' @param vector a character vector of genotypes
//' @param sep separator characters (for example "", "/", or "/|")
//' @return List of alleles, counts, frequencies, missing genotype count and rate.
//' @details Missing alleles (NA, '.', empty strings) are excluded from allele
//' counts. A genotype with either allele missing contributes to missing counts.
//' @export
// [[Rcpp::export]]
Rcpp::List stats_freq(Rcpp::CharacterVector& vector,std::string sep) {
    GWLD::GenotypeEncoder encoder(sep);
    encoder.analyze(STRING_PTR_RO(vector),vector.size());
    return encoder.frequencies(vector.size());
}

//' Convert diploid character genotypes to major-allele dosages
//' @param vector a character vector of genotypes
//' @param sep separator characters (for example "", "/", or "/|")
//' @details Outputs 0 for two major alleles, 1 for one, 2 for none, and NA for
//' missing genotypes. Frequency ties use lexical allele order. Alleles are
//' matched as whole tokens, including multi-character alleles.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector recode(Rcpp::CharacterVector& vector,std::string sep) {
    Rcpp::IntegerVector output(vector.size());
    GWLD::GenotypeEncoder encoder(sep);
    encoder.analyze(STRING_PTR_RO(vector),vector.size(),output.begin());
    return output;
}

// [[Rcpp::export]]
arma::Mat<int> code_Mat(Rcpp::CharacterMatrix& Matrix,std::string sep) {
    const int rows=Matrix.nrow(),cols=Matrix.ncol();
    arma::Mat<int> output(rows,cols);
    if (!rows || !cols) return output;
    const SEXP* strings=STRING_PTR_RO(Matrix);
    GWLD::GenotypeEncoder encoder(sep);
    for (int column=0;column<cols;++column) {
        encoder.analyze(strings+static_cast<R_xlen_t>(column)*rows,rows,output.colptr(column));
        if (column%256==0) Rcpp::checkUserInterrupt();
    }
    return output;
}
