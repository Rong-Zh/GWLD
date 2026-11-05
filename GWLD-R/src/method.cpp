#define ARMA_64BIT_WORD 1
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "../inst/include/gwld/engine.hpp"

namespace {
double scalar(arma::Col<int>& a,arma::Col<int>& b,GWLD::Method method) {
    if (a.n_elem!=b.n_elem) Rcpp::stop("g1 and g2 have different length");
    const auto sa=GWLD::summarize(a.memptr(),a.n_elem),sb=GWLD::summarize(b.memptr(),b.n_elem);
    const double value=GWLD::pair(a.memptr(),b.memptr(),a.n_elem,sa,sb,method);
    return std::isnan(value) ? NA_REAL : value;
}
void threads_valid(int threads) {
    if (threads<1) Rcpp::stop("Thread count must be positive");
    #ifndef _OPENMP
    if (threads>1) Rcpp::stop("This build has no OpenMP support; use one thread");
    #endif
}
arma::Mat<double> matrix(arma::Mat<int>& g,const std::string& method,int threads) {
    threads_valid(threads);
    const auto selected=GWLD::parse_method(method);
    return GWLD::Engine(g).matrix(selected,threads,NA_REAL);
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
Rcpp::CharacterMatrix collect(const std::vector<std::vector<Result>>& rows,const std::vector<GWLD::Site>& sites,bool decay) {
    std::size_t size=0;
    for (const auto& row:rows) size+=row.size();
    if (size>static_cast<std::size_t>(std::numeric_limits<int>::max())) Rcpp::stop("Too many output pairs for an R matrix");
    Rcpp::CharacterMatrix result(size,decay ? 4 : 7);
    std::size_t index=0;
    for (const auto& row:rows) for (const auto& r:row) {
        const auto& a=sites[r.i]; const auto& b=sites[r.j];
        if (decay) {
            result(index,0)=a.id; result(index,1)=b.id;
            result(index,2)=std::to_string(std::llabs(b.pos-a.pos));
        } else {
            result(index,0)=a.chr; result(index,1)=std::to_string(a.pos); result(index,2)=a.id;
            result(index,3)=b.chr; result(index,4)=std::to_string(b.pos); result(index,5)=b.id;
        }
        result(index,decay ? 3 : 6)=std::isnan(r.value) ? "NA" : std::to_string(r.value);
        ++index;
    }
    Rcpp::colnames(result)=decay ? Rcpp::CharacterVector::create("ID_1","ID_2","Dist","Value")
        : Rcpp::CharacterVector::create("CHR_1","POS_1","ID_1","CHR_2","POS_2","ID_2","Value");
    return result;
}
template<class Evaluate>
Rcpp::CharacterMatrix compute_pairwise(Rcpp::CharacterMatrix& info,std::size_t n,Evaluate evaluate,
                                  int threads,bool decay,int maxdist,double threshold) {
    threads_valid(threads);
    if (maxdist<0) Rcpp::stop("maxdist must be nonnegative");
    if (!std::isfinite(threshold)) Rcpp::stop("threshold must be finite");
    const auto sites=metadata(info,n);
    const GWLD::Window window(sites);
    std::vector<std::vector<Result>> rows(n);
    // Convert R strings before entering OpenMP; only C++ data in workers.
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for (std::int64_t i=0; i<static_cast<std::int64_t>(n); ++i) {
        auto emit=[&](std::size_t j) {
            if (!decay && sites[i].chr==sites[j].chr) return;
            const double value=evaluate(i,j);
            if (decay || value>threshold) rows[i].push_back({static_cast<std::size_t>(i),j,value});
        };
        if (decay) for (auto j:window.partners(i,sites[i],static_cast<std::int64_t>(maxdist)*1000)) emit(j);
        else for (std::size_t j=i+1; j<n; ++j) emit(j);
    }
    Rcpp::checkUserInterrupt();
    return collect(rows,sites,decay);
}
void square(const arma::Mat<double>& m) {
    if (m.n_rows!=m.n_cols) Rcpp::stop("Expected a square upper-triangular matrix");
}
} // namespace

// [[Rcpp::export]]
double RMIC(arma::Col<int>& g1,arma::Col<int>& g2) { return scalar(g1,g2,GWLD::Method::RMI); }
// [[Rcpp::export]]
arma::Mat<double> RMIC_Mat(arma::Mat<int>& geno012,int cores=1) { return matrix(geno012,"RMI",cores); }
// [[Rcpp::export]]
double MIC(arma::Col<int>& g1,arma::Col<int>& g2) { return scalar(g1,g2,GWLD::Method::MI); }
// [[Rcpp::export]]
arma::Mat<double> MIC_Mat(arma::Mat<int>& geno012,int cores=1) { return matrix(geno012,"MI",cores); }
// [[Rcpp::export]]
double LDC(arma::Col<int>& g1,arma::Col<int>& g2,std::string method="r^2") {
    const auto selected=GWLD::parse_method(method);
    if (selected==GWLD::Method::MI || selected==GWLD::Method::RMI) Rcpp::stop("LD requires r^2, r2, D or D'");
    return scalar(g1,g2,selected);
}
// [[Rcpp::export]]
arma::Mat<double> LDC_Mat(arma::Mat<int>& geno012,std::string method,int cores=1) {
    const auto selected=GWLD::parse_method(method);
    if (selected==GWLD::Method::MI || selected==GWLD::Method::RMI) Rcpp::stop("LD requires r^2, r2, D or D'");
    return matrix(geno012,method,cores);
}
// [[Rcpp::export]]
arma::Mat<double> gwld(arma::Mat<int>& geno012,std::string method="RMI",int thread=1) { return matrix(geno012,method,thread); }
// [[Rcpp::export]]
Rcpp::CharacterMatrix decay(arma::Mat<double>& m,Rcpp::CharacterMatrix& snpinfo,int maxdist=300) {
    square(m);
    return compute_pairwise(snpinfo,m.n_cols,[&](std::size_t i,std::size_t j){return m(i,j);},1,true,maxdist,0);
}
// [[Rcpp::export]]
Rcpp::CharacterMatrix calc_decay(arma::Mat<int>& geno012,Rcpp::CharacterMatrix& snpinfo,std::string method="RMI",int thread=1,int maxdist=300) {
    const auto selected=GWLD::parse_method(method); const GWLD::Engine engine(geno012);
    return compute_pairwise(snpinfo,geno012.n_cols,[&](std::size_t i,std::size_t j){return engine.at(i,j,selected);},thread,true,maxdist,0);
}
// [[Rcpp::export]]
Rcpp::CharacterMatrix calc_circos(arma::Mat<int>& geno012,Rcpp::CharacterMatrix& snpinfo,std::string method="RMI",int thread=1,double threshold=0.2) {
    const auto selected=GWLD::parse_method(method); const GWLD::Engine engine(geno012);
    return compute_pairwise(snpinfo,geno012.n_cols,[&](std::size_t i,std::size_t j){return engine.at(i,j,selected);},thread,false,0,threshold);
}
// [[Rcpp::export]]
Rcpp::CharacterMatrix circos(arma::Mat<double>& m,Rcpp::CharacterMatrix& snpinfo,double threshold=0.2) {
    square(m);
    return compute_pairwise(snpinfo,m.n_cols,[&](std::size_t i,std::size_t j){return m(i,j);},1,false,0,threshold);
}
