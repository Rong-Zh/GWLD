#define ARMA_64BIT_WORD 1
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "../inst/include/gwld/engine.hpp"
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>

namespace {
Rcpp::IntegerVector checked_genotypes(SEXP input,bool matrix_input) {
    if (matrix_input && !Rf_isMatrix(input)) Rcpp::stop("Expected a numeric genotype matrix");
    if (TYPEOF(input)!=INTSXP && TYPEOF(input)!=REALSXP)
        Rcpp::stop("Genotypes must be numeric 0, 1, 2 or missing (-1/NA)");
    if (TYPEOF(input)==INTSXP) return Rcpp::IntegerVector(input); // Engine validates integer codes without copying.
    Rcpp::NumericVector source(input);
    Rcpp::IntegerVector result(source.size());
    for (R_xlen_t i=0;i<source.size();++i) {
        const double value=source[i];
        if (std::isnan(value)) result[i]=NA_INTEGER;
        else if (value==0 || value==1 || value==2 || value==-1) result[i]=static_cast<int>(value);
        else Rcpp::stop("Genotypes must be exact 0, 1, 2 or missing (-1/NA); invalid value at position %lld",static_cast<long long>(i+1));
        if (i%65536==0) Rcpp::checkUserInterrupt();
    }
    return result;
}
arma::Mat<int> matrix_view(Rcpp::IntegerVector& storage,SEXP original) {
    Rcpp::IntegerVector dims=Rf_getAttrib(original,R_DimSymbol);
    return arma::Mat<int>(storage.begin(),dims[0],dims[1],false,true);
}
void checkpoint() { Rcpp::checkUserInterrupt(); }
double scalar(arma::Col<int>& a,arma::Col<int>& b,GWLD::Method method) {
    if (a.n_elem!=b.n_elem) Rcpp::stop("g1 and g2 have different length");
    const auto sa=GWLD::summarize(a.memptr(),a.n_elem),sb=GWLD::summarize(b.memptr(),b.n_elem);
    const double value=GWLD::pair(a.memptr(),b.memptr(),a.n_elem,sa,sb,method);
    return std::isnan(value) ? NA_REAL : value;
}
void threads_valid(int threads) {
    if (threads<1) Rcpp::stop("Thread count must be positive");
    #if !defined(_OPENMP) || defined(GWLD_NO_OPENMP)
    if (threads>1) Rcpp::stop("This build has no OpenMP support; use one thread");
    #endif
}
arma::Mat<double> matrix(arma::Mat<int>& g,const std::string& method,int threads) {
    threads_valid(threads);
    const auto selected=GWLD::parse_method(method);
    return GWLD::Engine(g,selected).matrix(selected,threads,NA_REAL,checkpoint);
}
std::vector<GWLD::Site> metadata(Rcpp::CharacterMatrix& info,std::size_t n) {
    if (info.ncol()<3 || static_cast<std::size_t>(info.nrow())!=n)
        Rcpp::stop("snpinfo must have one row per variant and at least three columns (chr, pos, id)");
    std::vector<GWLD::Site> sites; sites.reserve(n);
    for (std::size_t i=0; i<n; ++i) {
        for (int j=0; j<3; ++j) if (Rcpp::CharacterVector::is_na(info(i,j))) Rcpp::stop("snpinfo contains NA");
        sites.push_back({Rcpp::as<std::string>(info(i,0)),Rcpp::as<std::string>(info(i,2)),GWLD::position(Rcpp::as<std::string>(info(i,1)))});
    }
    return sites;
}
struct Result { std::size_t i,j; double value; };
std::string numeric_text(double value) {
    if (std::isnan(value)) return "NA";
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text<<std::setprecision(std::numeric_limits<double>::max_digits10)<<value;
    return text.str();
}
Rcpp::CharacterMatrix collect(const std::vector<Result>& rows,const std::vector<GWLD::Site>& sites,bool decay) {
    Rcpp::CharacterMatrix result(rows.size(),decay ? 4 : 7);
    std::size_t index=0;
    for (const auto& r:rows) {
        const auto& a=sites[r.i]; const auto& b=sites[r.j];
        if (decay) {
            result(index,0)=a.id; result(index,1)=b.id;
            result(index,2)=std::to_string(std::llabs(b.pos-a.pos));
        } else {
            result(index,0)=a.chr; result(index,1)=std::to_string(a.pos); result(index,2)=a.id;
            result(index,3)=b.chr; result(index,4)=std::to_string(b.pos); result(index,5)=b.id;
        }
        result(index,decay ? 3 : 6)=numeric_text(r.value);
        ++index;
        if (index%4096==0) Rcpp::checkUserInterrupt();
    }
    Rcpp::colnames(result)=decay ? Rcpp::CharacterVector::create("ID_1","ID_2","Dist","Value")
        : Rcpp::CharacterVector::create("CHR_1","POS_1","ID_1","CHR_2","POS_2","ID_2","Value");
    return result;
}
template<class Evaluate>
Rcpp::RObject compute_pairwise(Rcpp::CharacterMatrix& info,std::size_t n,Evaluate evaluate,
                               int threads,bool decay,int maxdist,double threshold,const std::string& outputfile) {
    threads_valid(threads);
    if (maxdist<0) Rcpp::stop("maxdist must be nonnegative");
    if (!std::isfinite(threshold)) Rcpp::stop("threshold must be finite");
    const auto sites=metadata(info,n);
    const GWLD::Window window(sites);
    std::ofstream stream;
    if (!outputfile.empty()) {
        stream.open(outputfile);
        if (!stream) Rcpp::stop("Cannot create pair output file: %s",outputfile.c_str());
        stream.exceptions(std::ios::failbit|std::ios::badbit);
        stream.imbue(std::locale::classic());
        stream<<std::setprecision(std::numeric_limits<double>::max_digits10);
        stream<<(decay ? "ID_1\tID_2\tDist\tValue\n" : "CHR_1\tPOS_1\tID_1\tCHR_2\tPOS_2\tID_2\tValue\n");
    }
    std::vector<Result> collected;
    const std::size_t batch=std::max<std::size_t>(1,std::min<std::size_t>(32,65536/std::max<std::size_t>(1,n)));
    for (std::size_t start=0;start<n;start+=batch) {
        Rcpp::checkUserInterrupt();
        const auto count=std::min(batch,n-start);
        std::vector<std::vector<Result>> rows(count);
        // R APIs and file writes stay on the calling thread, outside OpenMP.
        #if defined(_OPENMP) && !defined(GWLD_NO_OPENMP)
        #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        #endif
        for (std::int64_t offset=0;offset<static_cast<std::int64_t>(count);++offset) {
            const auto i=start+offset;
            auto emit=[&](std::size_t j) {
                if (!decay && sites[i].chr==sites[j].chr) return;
                const double value=evaluate(i,j);
                if (decay || value>threshold) rows[offset].push_back({i,j,value});
            };
            if (decay) for (auto j:window.partners(i,sites[i],static_cast<std::int64_t>(maxdist)*1000)) emit(j);
            else for (std::size_t j=i+1;j<n;++j) emit(j);
        }
        for (const auto& row:rows) {
            if (outputfile.empty()) {
                if (row.size()>static_cast<std::size_t>(std::numeric_limits<int>::max())-collected.size())
                    Rcpp::stop("Too many output pairs for an R matrix; use file= to stream results");
                collected.insert(collected.end(),row.begin(),row.end());
            } else for (const auto& pair:row) {
                const auto& a=sites[pair.i]; const auto& b=sites[pair.j];
                if (decay) stream<<a.id<<'\t'<<b.id<<'\t'<<std::llabs(b.pos-a.pos);
                else stream<<a.chr<<'\t'<<a.pos<<'\t'<<a.id<<'\t'<<b.chr<<'\t'<<b.pos<<'\t'<<b.id;
                stream<<'\t';
                if (std::isnan(pair.value)) stream<<"NA"; else stream<<pair.value;
                stream<<'\n';
            }
        }
    }
    Rcpp::checkUserInterrupt();
    if (!outputfile.empty()) { stream.close(); return Rcpp::wrap(outputfile); }
    return collect(collected,sites,decay);
}
void square(const arma::Mat<double>& m) {
    if (m.n_rows!=m.n_cols) Rcpp::stop("Expected a square upper-triangular matrix");
}
} // namespace

// [[Rcpp::export]]
double RMIC(SEXP g1,SEXP g2) {
    auto aa=checked_genotypes(g1,false),bb=checked_genotypes(g2,false);
    arma::Col<int> a(aa.begin(),aa.size(),false,true),b(bb.begin(),bb.size(),false,true);
    return scalar(a,b,GWLD::Method::RMI);
}
// [[Rcpp::export]]
arma::Mat<double> RMIC_Mat(SEXP geno012,int cores=1) {
    auto storage=checked_genotypes(geno012,true); auto g=matrix_view(storage,geno012);
    return matrix(g,"RMI",cores);
}
// [[Rcpp::export]]
double MIC(SEXP g1,SEXP g2) {
    auto aa=checked_genotypes(g1,false),bb=checked_genotypes(g2,false);
    arma::Col<int> a(aa.begin(),aa.size(),false,true),b(bb.begin(),bb.size(),false,true);
    return scalar(a,b,GWLD::Method::MI);
}
// [[Rcpp::export]]
arma::Mat<double> MIC_Mat(SEXP geno012,int cores=1) {
    auto storage=checked_genotypes(geno012,true); auto g=matrix_view(storage,geno012);
    return matrix(g,"MI",cores);
}
// [[Rcpp::export]]
double LDC(SEXP g1,SEXP g2,std::string method="r^2") {
    const auto selected=GWLD::parse_method(method);
    if (selected==GWLD::Method::MI || selected==GWLD::Method::RMI) Rcpp::stop("LD requires r^2, r2, D or D'");
    auto aa=checked_genotypes(g1,false),bb=checked_genotypes(g2,false);
    arma::Col<int> a(aa.begin(),aa.size(),false,true),b(bb.begin(),bb.size(),false,true);
    return scalar(a,b,selected);
}
// [[Rcpp::export]]
arma::Mat<double> LDC_Mat(SEXP geno012,std::string method,int cores=1) {
    const auto selected=GWLD::parse_method(method);
    if (selected==GWLD::Method::MI || selected==GWLD::Method::RMI) Rcpp::stop("LD requires r^2, r2, D or D'");
    auto storage=checked_genotypes(geno012,true); auto g=matrix_view(storage,geno012);
    return matrix(g,method,cores);
}
// [[Rcpp::export]]
arma::Mat<double> gwld(SEXP geno012,std::string method="RMI",int thread=1) {
    auto storage=checked_genotypes(geno012,true); auto g=matrix_view(storage,geno012);
    return matrix(g,method,thread);
}
// [[Rcpp::export]]
Rcpp::RObject decay(arma::Mat<double>& m,Rcpp::CharacterMatrix& snpinfo,int maxdist=300,std::string outputfile="") {
    square(m);
    return compute_pairwise(snpinfo,m.n_cols,[&](std::size_t i,std::size_t j){return m(i,j);},1,true,maxdist,0,outputfile);
}
// [[Rcpp::export]]
Rcpp::RObject calc_decay(SEXP geno012,Rcpp::CharacterMatrix& snpinfo,std::string method="RMI",int thread=1,int maxdist=300,std::string outputfile="") {
    threads_valid(thread);
    auto storage=checked_genotypes(geno012,true); auto g=matrix_view(storage,geno012);
    const auto selected=GWLD::parse_method(method); const GWLD::Engine engine(g,selected);
    return compute_pairwise(snpinfo,g.n_cols,[&](std::size_t i,std::size_t j){return engine.at(i,j,selected);},thread,true,maxdist,0,outputfile);
}
// [[Rcpp::export]]
Rcpp::RObject calc_circos(SEXP geno012,Rcpp::CharacterMatrix& snpinfo,std::string method="RMI",int thread=1,double threshold=0.2,std::string outputfile="") {
    threads_valid(thread);
    auto storage=checked_genotypes(geno012,true); auto g=matrix_view(storage,geno012);
    const auto selected=GWLD::parse_method(method); const GWLD::Engine engine(g,selected);
    return compute_pairwise(snpinfo,g.n_cols,[&](std::size_t i,std::size_t j){return engine.at(i,j,selected);},thread,false,0,threshold,outputfile);
}
// [[Rcpp::export]]
Rcpp::RObject circos(arma::Mat<double>& m,Rcpp::CharacterMatrix& snpinfo,double threshold=0.2,std::string outputfile="") {
    square(m);
    return compute_pairwise(snpinfo,m.n_cols,[&](std::size_t i,std::size_t j){return m(i,j);},1,false,0,threshold,outputfile);
}
