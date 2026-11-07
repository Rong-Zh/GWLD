#include "utils.h"
#include "gwld/core.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <cstdint>
#include <fstream>
#include <limits>
#include <locale>
#include <stdexcept>
#include <unordered_map>

// Keep C++11 compatibility; use allocation-free numeric formatting with recent libraries.
#if __cplusplus >= 201703L && ((defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE >= 11) || \
    (defined(_LIBCPP_VERSION) && _LIBCPP_VERSION >= 170000 && \
     (!defined(__APPLE__) || __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ >= 140000)))
#include <charconv>
#define GWLD_UTILS_CHARCONV 1
#endif

namespace {
class Writer {
    std::vector<char> file_buffer;
    std::ofstream out;
    std::string buffer;
public:
    explicit Writer(const std::string& path) : file_buffer(65536) {
        out.rdbuf()->pubsetbuf(file_buffer.data(),file_buffer.size());
        out.open(path);
        if (!out) throw std::runtime_error("Cannot create output: "+path);
        out.exceptions(std::ios::failbit|std::ios::badbit);
        out.imbue(std::locale::classic());
        buffer.reserve(1u<<20);
    }
    void flush() {
        if (!buffer.empty()) { out.write(buffer.data(),static_cast<std::streamsize>(buffer.size())); buffer.clear(); }
    }
    void value(const std::string& text) {
        if (buffer.size()+text.size()>(1u<<20)) flush();
        if (text.size()>(1u<<20)) out.write(text.data(),static_cast<std::streamsize>(text.size()));
        else buffer.append(text);
    }
    void value(double value) {
        #if defined(GWLD_UTILS_CHARCONV)
        if (std::isnan(value)) {
            static const std::array<std::string,2> names=[]() {
                std::array<std::string,2> result;
                for (unsigned i=0;i<2;++i) {
                    std::ostringstream stream; stream.imbue(std::locale::classic());
                    stream<<(i ? -std::numeric_limits<double>::quiet_NaN() : std::numeric_limits<double>::quiet_NaN());
                    result[i]=stream.str();
                }
                return result;
            }();
            buffer.append(names[std::signbit(value) ? 1:0]); return;
        }
        char text[64];
        const auto result=std::to_chars(text,text+sizeof(text),value,std::chars_format::general,6);
        if (result.ec!=std::errc()) throw std::runtime_error("Numeric formatting failed");
        buffer.append(text,result.ptr);
        #else
        // Reuse the stream's formatter instead of constructing a stream per number.
        flush(); out<<value;
        #endif
    }
    void character(char ch) { buffer.push_back(ch); if (buffer.size()>=(1u<<20)) flush(); }
    void row_end() { character('\n'); if (buffer.size()>=(1u<<20)) flush(); }
    void close() { flush(); out.close(); }
};
std::array<double,3> frequencies(const int* values,arma::uword size) {
    const auto summary=GWLD::summarize(values,size);
    return {{summary.p,1-summary.p,size ? 1-summary.total/static_cast<double>(size) : std::numeric_limits<double>::quiet_NaN()}};
}
}

std::vector<char> is_any_of(const std::string& str) { return {str.begin(),str.end()}; }
void split(std::vector<std::string>& result,const std::string& str,const std::vector<char>& delimiters) {
    // Preserve the old by-value behavior when the source is an output element.
    if (!result.empty()) {
        const auto address=reinterpret_cast<std::uintptr_t>(&str);
        if (address>=reinterpret_cast<std::uintptr_t>(result.data()) &&
            address<reinterpret_cast<std::uintptr_t>(result.data()+result.size())) {
            const std::string saved=str;
            split(result,saved,delimiters); return;
        }
    }
    bool mask[256]={};
    for (unsigned char ch:delimiters) mask[ch]=true;
    std::size_t start=0;
    for (std::size_t i=0;i<str.size();++i) if (mask[static_cast<unsigned char>(str[i])]) {
        result.emplace_back(str,start,i-start); start=i+1;
    }
    if (start<str.size()) result.emplace_back(str,start,str.size()-start);
}
std::vector<std::vector<std::string>> readfile(const std::string& filename,char sep) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open input: "+filename);
    std::vector<std::vector<std::string>> result;
    std::string line;
    while (std::getline(in,line)) {
        std::vector<std::string> row;
        const auto fields=static_cast<std::size_t>(std::count(line.begin(),line.end(),sep))+
            static_cast<std::size_t>(!line.empty() && line.back()!=sep);
        row.reserve(fields);
        std::size_t start=0;
        while (start<line.size()) {
            const auto end=line.find(sep,start);
            row.emplace_back(line,start,end==std::string::npos ? line.size()-start : end-start);
            if (end==std::string::npos) break;
            start=end+1;
        }
        result.push_back(std::move(row));
    }
    if (!in.eof()) throw std::runtime_error("Cannot read input: "+filename);
    return result;
}
template<class T> void writefile(const T& data,const std::string& filename,char sep) {
    Writer out(filename);
    for (const auto& row:data) {
        for (const auto& value:row) { out.value(value); out.character(sep); }
        out.row_end();
    }
    out.close();
}
template void writefile(const std::vector<std::vector<std::string>>&,const std::string&,char);
template void writefile(const std::vector<std::vector<double>>&,const std::string&,char);
void writearma(const arma::Mat<double>& data,const std::string& filename,char sep) {
    Writer out(filename);
    for (arma::uword row=0;row<data.n_rows;++row) {
        for (arma::uword col=0;col<data.n_cols;++col) { out.value(data(row,col)); out.character(sep); }
        out.row_end();
    }
    out.close();
}
bool cmp(const std::pair<std::string,int>& a,const std::pair<std::string,int>& b) {
    return a.second!=b.second ? a.second>b.second : a.first<b.first;
}
std::vector<std::pair<std::string,int>> counts_allele(const std::vector<std::string>& alleles) {
    if (alleles.empty()) return {};
    std::unordered_map<std::string,std::size_t> counts;
    counts.reserve(8);
    for (const auto& allele:alleles) ++counts[allele];
    std::vector<std::pair<std::string,int>> result; result.reserve(counts.size());
    for (const auto& item:counts) {
        if (item.second>static_cast<std::size_t>(std::numeric_limits<int>::max())) throw std::overflow_error("Allele count exceeds integer range");
        result.emplace_back(item.first,static_cast<int>(item.second));
    }
    std::sort(result.begin(),result.end(),cmp);
    return result;
}
std::vector<double> allele_freq(const arma::Col<int>& genotype) {
    const auto result=frequencies(genotype.memptr(),genotype.n_elem);
    return std::vector<double>(result.begin(),result.end());
}
arma::Mat<double> table(const arma::Col<int>& a,const arma::Col<int>& b,int rows,int cols) {
    if (a.n_elem!=b.n_elem || rows<0 || cols<0) throw std::invalid_argument("Invalid contingency-table dimensions");
    arma::Mat<double> result(rows,cols,arma::fill::zeros);
    std::array<std::uint64_t,9> small{};
    const bool fixed=rows==3 && cols==3;
    for (arma::uword i=0;i<a.n_elem;++i) {
        const int x=a[i],y=b[i];
        if ((x<0 && x!=-1 && x!=std::numeric_limits<int>::min()) ||
            (y<0 && y!=-1 && y!=std::numeric_limits<int>::min())) throw std::invalid_argument("Invalid genotype");
        if (x>=rows || y>=cols) throw std::invalid_argument("Genotype exceeds contingency-table dimensions");
        if (x<0 || y<0) continue;
        if (fixed) ++small[x+3*y]; else ++result.memptr()[x+static_cast<arma::uword>(rows)*y];
    }
    if (fixed) for (unsigned i=0;i<9;++i) result[i]=static_cast<double>(small[i]);
    return result;
}
std::vector<std::vector<double>> stats_freq(const arma::Mat<int>& genotypes) {
    std::vector<std::vector<double>> result; result.reserve(genotypes.n_cols);
    for (arma::uword j=0;j<genotypes.n_cols;++j) {
        const auto values=frequencies(genotypes.colptr(j),genotypes.n_rows);
        result.emplace_back(values.begin(),values.end());
    }
    return result;
}
std::vector<std::vector<std::string>> getsnpinfo(const std::vector<std::vector<std::string>>& data,const std::vector<int>& columns) {
    std::vector<std::vector<std::string>> result; result.reserve(data.size());
    for (const auto& row:data) {
        std::vector<std::string> selected; selected.reserve(columns.size());
        for (int column:columns) {
            if (column<0 || static_cast<std::size_t>(column)>=row.size()) throw std::invalid_argument("SNP metadata column is out of range");
            selected.push_back(row[column]);
        }
        result.push_back(std::move(selected));
    }
    return result;
}
std::vector<std::vector<std::string>> getchrinfo(const std::vector<std::vector<std::string>>& sites) {
    std::vector<std::vector<std::string>> result;
    std::size_t start=0;
    while (start<sites.size()) {
        if (sites[start].empty()) throw std::invalid_argument("Missing chromosome metadata");
        std::size_t end=start+1;
        while (end<sites.size() && !sites[end].empty() && sites[end][0]==sites[start][0]) ++end;
        result.push_back({sites[start][0],std::to_string(start),std::to_string(end-1),std::to_string(end-start)});
        start=end;
    }
    return result;
}
