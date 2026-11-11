#pragma once
#include "packed.hpp"

// Include after armadillo or RcppArmadillo. Workers touch only C++ memory.
namespace GWLD {
class Engine {
    const arma::Mat<int>& genotypes;
    bool has_constant=false;
    std::vector<Summary> summaries;
    PackedGenotypes packed;
    std::vector<double> log_factorials,count_logs,half_logs,mi_marginals;
    struct RmiMarginal { double log_sum=0,shape=0,gamma=0; };
    struct RmiLocus {
        bool complete=false;
        unsigned categories=0;
        double factorial_sum=0;
        std::array<RmiMarginal,3> marginal{};
    };
    std::vector<RmiLocus> rmi_loci;
    double rmi_normalizers[3][3]{};
    double factorial(std::size_t count) const {
        return log_factorials.empty() ? log_gamma(static_cast<double>(count)+1) : log_factorials[count];
    }
    void prepare_rmi() {
        if (genotypes.n_cols<4 || !genotypes.n_rows) return;
        const double n=static_cast<double>(genotypes.n_rows);
        rmi_loci.resize(summaries.size());
        for (unsigned r=1;r<=3;++r) for (unsigned s=1;s<=3;++s)
            rmi_normalizers[r-1][s-1]=(r-1)*(s-1)*std::log(n+0.5*r*s);
        for (std::size_t j=0;j<summaries.size();++j) {
            const auto& summary=summaries[j]; auto& locus=rmi_loci[j];
            if (summary.total!=n) continue;
            locus.complete=true;
            for (unsigned k=0;k<3;++k) if (summary.present[k]) {
                ++locus.categories;
                locus.factorial_sum+=factorial(static_cast<std::size_t>(summary.counts[k]));
            }
            const double r=locus.categories;
            for (unsigned s=1;s<=3;++s) {
                auto& marginal=locus.marginal[s-1];
                const double w=n/(n+0.5*r*s); double square_sum=0;
                for (unsigned k=0;k<3;++k) if (summary.present[k]) {
                    const double x=(1-w)/r+w*summary.counts[k]/n;
                    square_sum+=x*x; marginal.log_sum+=std::log(x);
                }
                const double shape=(s+1)/(s*square_sum)-1.0/s;
                marginal.shape=shape;
                marginal.gamma=log_gamma(shape*s)-s*log_gamma(shape)-r*log_gamma(s);
            }
        }
    }
    double complete_rmi(const Counts& counts,std::size_t i,std::size_t j) const {
        const auto& a=rmi_loci[i]; const auto& b=rmi_loci[j];
        const auto& x=a.marginal[b.categories-1]; const auto& y=b.marginal[a.categories-1];
        double value=factorial(genotypes.n_rows);
        for (double cell:counts.cells) value+=factorial(static_cast<std::size_t>(cell));
        value-=a.factorial_sum+b.factorial_sum;
        const double omega=rmi_normalizers[a.categories-1][b.categories-1]
            +0.5*(a.categories+x.shape-2)*y.log_sum+0.5*(b.categories+y.shape-2)*x.log_sum
            +0.5*(x.gamma+y.gamma);
        return (value-omega)/(counts.n*std::log(2.0));
    }
public:
    Engine(const arma::Mat<int>& g,Method method=Method::RMI) : genotypes(g) {
        summaries.reserve(g.n_cols);
        if (g.n_rows>=128 && g.n_cols>=4) {
            packed.build(g);
            // Reuse packed diagonals instead of scanning and validating every
            // integer genotype a second time just to obtain locus marginals.
            for (arma::uword i=0;i<g.n_cols;++i) {
                const auto counts=packed.count(i,i); Summary summary;
                summary.total=counts.n;
                for (unsigned k=0;k<3;++k) {
                    summary.counts[k]=counts.rows[k]; summary.present[k]=counts.rows[k]>0;
                    if (counts.rows[k]==g.n_rows) summary.constant=static_cast<signed char>(k);
                }
                if (counts.n) summary.p=(counts.rows[0]+0.5*counts.rows[1])/counts.n;
                summaries.push_back(summary);
            }
        } else for (arma::uword i=0;i<g.n_cols;++i)
            summaries.push_back(summarize(g.n_rows ? g.colptr(i) : nullptr,g.n_rows));
        has_constant=std::any_of(summaries.begin(),summaries.end(),[](const Summary& s) { return s.constant>=0; });
        // Cache only when repeated pair evaluations amortize construction, and
        // bound the cache to 8 MiB even for very large sample cohorts.
        if (method==Method::RMI && g.n_rows>0 && g.n_rows<(1u<<20) &&
            static_cast<double>(g.n_cols)*(g.n_cols ? g.n_cols-1 : 0)*8>=g.n_rows+1) {
            log_factorials.resize(g.n_rows+1);
            for (std::size_t i=0;i<log_factorials.size();++i)
                log_factorials[i]=log_gamma(static_cast<double>(i)+1);
        }
        if (method==Method::MI && g.n_rows>0 && g.n_rows<(1u<<20) &&
            static_cast<double>(g.n_cols)*(g.n_cols ? g.n_cols-1 : 0)*8>=g.n_rows+1) {
            count_logs.resize(g.n_rows+1);
            for (std::size_t i=1;i<count_logs.size();++i) count_logs[i]=i*std::log2(static_cast<double>(i));
            mi_marginals.resize(summaries.size());
            for (std::size_t j=0;j<summaries.size();++j) if (summaries[j].total==g.n_rows)
                for (unsigned k=0;k<3;++k) mi_marginals[j]+=count_logs[static_cast<std::size_t>(summaries[j].counts[k])];
        }
        if (method==Method::RMI) {
            prepare_rmi();
            // Integer/half-integer log arguments recur across missing-data pairs.
            // Use the original log function to tabulate them, not an approximation.
            // Allocate only for partial columns and cap this extra table at 8 MiB.
            if (!log_factorials.empty() && g.n_rows<((1u<<19)-5) &&
                static_cast<double>(g.n_cols)*(g.n_cols-1)>=g.n_rows &&
                std::any_of(summaries.begin(),summaries.end(),[&](const Summary& s) { return s.total>0 && s.total<g.n_rows; })) {
                half_logs.resize(2*g.n_rows+10);
                for (std::size_t k=1;k<half_logs.size();++k) half_logs[k]=std::log(0.5*k);
            }
        }
    }
private:
    template<bool Constant>
    double at_impl(std::size_t i, std::size_t j, Method method) const {
        if (!genotypes.n_rows) return nan();
        const bool complete=summaries[i].total==genotypes.n_rows && summaries[j].total==genotypes.n_rows;
        const Counts counts=Constant && (summaries[i].constant>=0 || summaries[j].constant>=0) ?
            count_constant(summaries[i],summaries[j]) : packed.ready() ?
            (complete ? packed.count_complete(i,j,summaries[i],summaries[j]) : packed.count(i,j))
            : GWLD::count(genotypes.colptr(i),genotypes.colptr(j),genotypes.n_rows);
        if (method==Method::MI && !count_logs.empty()) {
            if (complete) {
                double value=count_logs[genotypes.n_rows];
                for (double cell:counts.cells) value+=count_logs[static_cast<std::size_t>(cell)];
                return std::max(0.0,(value-mi_marginals[i]-mi_marginals[j])/counts.n);
            }
            return information(counts,summaries[i],summaries[j],false,nullptr,count_logs.data());
        }
        if (method==Method::RMI && !rmi_loci.empty() && rmi_loci[i].complete && rmi_loci[j].complete)
            return complete_rmi(counts,i,j);
        if (method==Method::RMI && !half_logs.empty())
            return information(counts,summaries[i],summaries[j],true,log_factorials.data(),nullptr,half_logs.data());
        return evaluate_counts(counts,summaries[i],summaries[j],method,
                               log_factorials.empty() ? nullptr : log_factorials.data());
    }
    template<bool Constant>
    arma::Mat<double> matrix_impl(Method method, int threads, double missing,void (*checkpoint)()) const {
        if (threads<1) throw std::invalid_argument("Threads must be positive");
        arma::Mat<double> result(genotypes.n_cols,genotypes.n_cols,arma::fill::none);
        if (result.n_cols) std::fill(result.colptr(0),result.colptr(0)+result.n_rows,missing);
        const auto n=static_cast<std::int64_t>(genotypes.n_cols);
        // R callbacks run on the calling thread between bounded batches.
        // CLI callers keep one parallel region and pay no callback overhead.
        const std::int64_t batch=checkpoint ? 32 : std::max<std::int64_t>(1,n);
        for (std::int64_t start=1;start<n;start+=batch) {
            if (checkpoint) checkpoint();
            const auto end=std::min(n,start+batch);
            #if defined(_OPENMP) && !defined(GWLD_NO_OPENMP)
            #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
            #endif
            for (std::int64_t j=start; j<end; ++j) {
                // Each output element is written once: lower triangle stays missing.
                std::fill(result.colptr(j)+j,result.colptr(j)+result.n_rows,missing);
                for (std::int64_t i=0; i<j; ++i) {
                    const double value=at_impl<Constant>(i,j,method);
                    result(i,j)=std::isnan(value) ? missing : value;
                }
            }
        }
        if (checkpoint) checkpoint();
        return result;
    }
public:
    double at(std::size_t i,std::size_t j,Method method) const {
        return has_constant ? at_impl<true>(i,j,method) : at_impl<false>(i,j,method);
    }
    arma::Mat<double> matrix(Method method,int threads,double missing=nan(),void (*checkpoint)()=nullptr) const {
        // Select once, so ordinary matrices pay no per-pair constant checks.
        return has_constant ? matrix_impl<true>(method,threads,missing,checkpoint)
                            : matrix_impl<false>(method,threads,missing,checkpoint);
    }

};
} // namespace GWLD
