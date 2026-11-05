#include <armadillo>
#include "gwld/engine.hpp"
#include <chrono>
#include <iostream>
#include <random>

arma::Mat<double> methods_mat(arma::Mat<int>&,const std::string&,int);
int main() {
    arma::Mat<int> g(512,384);
    std::mt19937 random(2026);
    for (auto& value:g) value=static_cast<int>(random()%4)-1;
    std::cout<<"method\tv1_seconds\tv2_seconds\tspeedup\tmax_error\n";
    for (const std::string method:{"MI","RMI"}) {
        const auto start=std::chrono::steady_clock::now();
        const auto old=methods_mat(g,method,1);
        const auto middle=std::chrono::steady_clock::now();
        const auto fresh=GWLD::Engine(g).matrix(GWLD::parse_method(method),1);
        const auto end=std::chrono::steady_clock::now();
        double error=0;
        for (arma::uword j=1;j<g.n_cols;++j) for(arma::uword i=0;i<j;++i)
            error=std::max(error,std::abs(old(i,j)-fresh(i,j)));
        const double a=std::chrono::duration<double>(middle-start).count(),b=std::chrono::duration<double>(end-middle).count();
        std::cout<<method<<'\t'<<a<<'\t'<<b<<'\t'<<a/b<<'\t'<<error<<'\n';
        if (error>1e-10) return 1;
    }
}
