#pragma once
#include "core.hpp"

// Include after armadillo or RcppArmadillo. Workers touch only C++ memory.
namespace GWLD {
class Engine {
    const arma::Mat<int>& genotypes;
    std::vector<Summary> summaries;
public:
    Engine(const arma::Mat<int>& g) : genotypes(g) {
        summaries.reserve(g.n_cols);
        for (arma::uword i=0; i<g.n_cols; ++i)
            summaries.push_back(summarize(g.colptr(i),g.n_rows));
    }
    double at(std::size_t i, std::size_t j, Method method) const {
        return pair(genotypes.colptr(i),genotypes.colptr(j),genotypes.n_rows,
                    summaries[i],summaries[j],method);
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
