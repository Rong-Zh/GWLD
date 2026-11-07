#include <armadillo>
#include "gwld/engine.hpp"
#include <chrono>
#include <iostream>
#include <random>

// Scalar counting reference using the current summaries, formulas and solver.
// This isolates execution paths; it is not a frozen historical-version benchmark.
arma::Mat<double> scalar_matrix(const arma::Mat<int>& g,GWLD::Method method,int threads) {
    std::vector<GWLD::Summary> summaries;
    for (arma::uword j=0;j<g.n_cols;++j) summaries.push_back(GWLD::summarize(g.n_rows ? g.colptr(j) : nullptr,g.n_rows));
    arma::Mat<double> result(g.n_cols,g.n_cols); result.fill(GWLD::nan());
    #ifndef _OPENMP
    (void)threads;
    #endif
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for (std::int64_t j=1;j<static_cast<std::int64_t>(g.n_cols);++j)
        for (std::int64_t i=0;i<j;++i)
            result(i,j)=GWLD::pair(g.n_rows ? g.colptr(i) : nullptr,g.n_rows ? g.colptr(j) : nullptr,
                                   g.n_rows,summaries[i],summaries[j],method);
    return result;
}
int main(int argc,char* argv[]) {
    try {
        const std::size_t samples=argc>1 ? std::stoul(argv[1]) : 512,variants=argc>2 ? std::stoul(argv[2]) : 384;
        const int threads=argc>3 ? std::stoi(argv[3]) : 1,missing=argc>4 ? std::stoi(argv[4]) : 25;
        if (threads<1 || missing<0 || missing>100) throw std::invalid_argument("Invalid threads or missing percentage");
        arma::Mat<int> g(samples,variants); std::mt19937 random(2026);
        for (auto& value:g) value=static_cast<int>(random()%100)<missing ? -1 : static_cast<int>(random()%3);
        std::cout<<"method\tscalar_seconds\tpacked_seconds\tspeedup\tmax_error\n";
        for (const std::string name:{"MI","RMI","r2","D","D'"}) {
            const auto method=GWLD::parse_method(name);
            std::vector<double> old_times,new_times; double error=0;
            for (int repeat=0;repeat<3;++repeat) {
                arma::Mat<double> old,fresh;
                auto run=[&](bool scalar) {
                    const auto start=std::chrono::steady_clock::now();
                    if (scalar) old=scalar_matrix(g,method,threads);
                    else fresh=GWLD::Engine(g,method).matrix(method,threads);
                    const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                    (scalar ? old_times : new_times).push_back(seconds);
                };
                run(repeat%2==0); run(repeat%2!=0);
                for (arma::uword j=1;j<g.n_cols;++j) for (arma::uword i=0;i<j;++i) {
                    if (std::isnan(old(i,j)) && std::isnan(fresh(i,j))) continue;
                    if (!std::isfinite(old(i,j)) || !std::isfinite(fresh(i,j))) throw std::runtime_error("Missing-value mismatch");
                    error=std::max(error,std::abs(old(i,j)-fresh(i,j)));
                }
            }
            std::sort(old_times.begin(),old_times.end()); std::sort(new_times.begin(),new_times.end());
            std::cout<<name<<'\t'<<old_times[1]<<'\t'<<new_times[1]<<'\t'<<old_times[1]/new_times[1]<<'\t'<<error<<'\n';
            if (error>1e-12) return 1;
        }
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n'; return 1;}
}
