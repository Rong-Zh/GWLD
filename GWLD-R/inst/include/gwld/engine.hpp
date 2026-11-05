#pragma once
#include "packed.hpp"

// Include after armadillo or RcppArmadillo. Workers touch only C++ memory.
namespace GWLD {
class Engine {
    const arma::Mat<int>& genotypes;
    std::vector<Summary> summaries;
    PackedGenotypes packed;
    std::vector<double> log_factorials;
public:
    Engine(const arma::Mat<int>& g,Method method=Method::RMI) : genotypes(g) {
        summaries.reserve(g.n_cols);
        for (arma::uword i=0; i<g.n_cols; ++i)
            summaries.push_back(summarize(g.n_rows ? g.colptr(i) : nullptr,g.n_rows));
        if (g.n_rows>=128 && g.n_cols>=4) packed.build(g);
        // Cache only when repeated pair evaluations amortize construction, and
        // bound the cache to 8 MiB even for very large sample cohorts.
        if (method==Method::RMI && g.n_rows>0 && g.n_rows<(1u<<20) &&
            static_cast<double>(g.n_cols)*(g.n_cols ? g.n_cols-1 : 0)*8>=g.n_rows+1) {
            log_factorials.resize(g.n_rows+1);
            for (std::size_t i=0;i<log_factorials.size();++i)
                log_factorials[i]=log_gamma(static_cast<double>(i)+1);
        }
    }
    double at(std::size_t i, std::size_t j, Method method) const {
        if (!genotypes.n_rows) return nan();
        const Counts counts=packed.ready() ? packed.count(i,j)
            : GWLD::count(genotypes.colptr(i),genotypes.colptr(j),genotypes.n_rows);
        return evaluate_counts(counts,summaries[i],summaries[j],method,
                               log_factorials.empty() ? nullptr : log_factorials.data());
    }
    arma::Mat<double> matrix(Method method, int threads, double missing=nan()) const {
        if (threads<1) throw std::invalid_argument("Threads must be positive");
        arma::Mat<double> result(genotypes.n_cols,genotypes.n_cols);
        result.fill(missing);
        const auto n=static_cast<std::int64_t>(genotypes.n_cols);
        // Each column is contiguous; dynamic scheduling balances triangular work.
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        #endif
        for (std::int64_t j=1; j<n; ++j)
            for (std::int64_t i=0; i<j; ++i) {
                const double value=at(i,j,method);
                result(i,j)=std::isnan(value) ? missing : value;
            }
        return result;
    }
};
} // namespace GWLD
