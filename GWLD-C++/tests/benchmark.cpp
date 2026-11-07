#include <armadillo>
#include "gwld/engine.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>

arma::Mat<double> methods_mat(arma::Mat<int>&,const std::string&,int);

// Compare the frozen pre-rewrite 1.x implementation with the current engine.
// Input generation and result verification are outside the timed region.
int main(int argc,char** argv) {
    try {
        if (argc>5) throw std::invalid_argument("Usage: gwld_benchmark [samples=512] [variants=384] [threads=1] [missing_percent=25]");
        const std::size_t samples=argc>1 ? std::stoul(argv[1]) : 512;
        const std::size_t variants=argc>2 ? std::stoul(argv[2]) : 384;
        const int threads=argc>3 ? std::stoi(argv[3]) : 1;
        const int missing=argc>4 ? std::stoi(argv[4]) : 25;
        if (samples<2 || variants<2 || threads<1 || missing<0 || missing>=100)
            throw std::invalid_argument("Require at least two samples/variants, positive threads and 0..99% missing");
        arma::Mat<int> g(samples,variants);
        std::mt19937 random(2026);
        for (auto& value:g) value=random()%100<static_cast<unsigned>(missing) ? -1 : static_cast<int>(random()%3);
        std::cout<<"method\tsamples\tvariants\tmissing_percent\tthreads\tv1_seconds\tv2_seconds\tspeedup\tmax_abs_diff\tmissing_mismatch\n"<<std::setprecision(17);
        for (const std::string name:{"MI","RMI","r^2","D","D'"}) {
            const auto method=GWLD::parse_method(name);
            std::vector<double> old_times,new_times;
            double error=0; std::size_t missing_mismatch=0;
            for (int repeat=-1;repeat<3;++repeat) {
                arma::Mat<double> old,fresh;
                const auto run=[&](bool legacy) {
                    const auto start=std::chrono::steady_clock::now();
                    if (legacy) old=methods_mat(g,name,threads);
                    else fresh=GWLD::Engine(g,method).matrix(method,threads);
                    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                    if (repeat>=0) (legacy ? old_times : new_times).push_back(elapsed);
                };
                run(repeat%2==0); run(repeat%2!=0);
                std::size_t mismatch=0;
                for (arma::uword j=1;j<g.n_cols;++j) for (arma::uword i=0;i<j;++i) {
                    const double a=old(i,j),b=fresh(i,j);
                    if (std::isnan(a) && std::isnan(b)) continue;
                    if (!std::isfinite(a) || !std::isfinite(b)) { ++mismatch; continue; }
                    error=std::max(error,std::abs(a-b));
                }
                missing_mismatch=std::max(missing_mismatch,mismatch);
            }
            std::sort(old_times.begin(),old_times.end());
            std::sort(new_times.begin(),new_times.end());
            std::cout<<name<<'\t'<<samples<<'\t'<<variants<<'\t'<<missing<<'\t'<<threads
                     <<'\t'<<old_times[1]<<'\t'<<new_times[1]<<'\t'<<old_times[1]/new_times[1]
                     <<'\t'<<error<<'\t'<<missing_mismatch<<std::endl;
            // LD uses a corrected optimizer, so report its differences explicitly.
            // MI/RMI retain the historical formula and must remain equivalent.
            if ((method==GWLD::Method::MI || method==GWLD::Method::RMI) &&
                (missing_mismatch || error>1e-10)) return 1;
        }
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
