// Export raw observations and production results for an independent R oracle.
#include <armadillo>
#include "gwld/engine.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

arma::Mat<double> methods_mat(arma::Mat<int>&,const std::string&,int);

int main(int argc,char** argv) {
    try {
        if (argc!=3) throw std::invalid_argument("Usage: gwld_ld_validation_export cases.tsv genotypes.tsv");
        std::ofstream cases(argv[1]),observations(argv[2]);
        if (!cases || !observations) throw std::runtime_error("Cannot create validation files");
        cases.exceptions(std::ios::badbit|std::ios::failbit);
        observations.exceptions(std::ios::badbit|std::ios::failbit);
        cases<<std::setprecision(17);
        cases<<"case_id\tdataset\tselection\tvariant_i\tvariant_j\tp\tq\told_D\tnew_D\tnew_r2\tnew_Dprime";
        for (int k=0;k<9;++k) cases<<"\tc"<<k;
        cases<<'\n';
        observations<<"case_id\tsample\tg1\tg2\n";
        std::size_t id=0;
        const auto emit=[&](const arma::Mat<int>& g,std::size_t i,std::size_t j,double old_D,
                            const std::string& dataset,const std::string& selection) {
            ++id;
            const auto a=GWLD::summarize(g.colptr(i),g.n_rows),b=GWLD::summarize(g.colptr(j),g.n_rows);
            const auto c=GWLD::count(g.colptr(i),g.colptr(j),g.n_rows);
            const GWLD::Engine engine(g,GWLD::Method::D);
            cases<<id<<'\t'<<dataset<<'\t'<<selection<<'\t'<<i+1<<'\t'<<j+1
                 <<'\t'<<a.p<<'\t'<<b.p<<'\t'<<old_D<<'\t'<<engine.at(i,j,GWLD::Method::D)
                 <<'\t'<<engine.at(i,j,GWLD::Method::R2)<<'\t'<<engine.at(i,j,GWLD::Method::Dprime);
            for (double cell:c.cells) cases<<'\t'<<cell;
            cases<<'\n';
            for (arma::uword k=0;k<g.n_rows;++k) observations<<id<<'\t'<<k+1<<'\t'<<g(k,i)<<'\t'<<g(k,j)<<'\n';
        };
        for (auto config:{std::array<int,3>{{512,384,25}},std::array<int,3>{{512,768,0}},std::array<int,3>{{4096,128,1}}}) {
            arma::Mat<int> g(config[0],config[1]); std::mt19937 random(2026);
            for (auto& x:g) x=random()%100<static_cast<unsigned>(config[2]) ? -1 : static_cast<int>(random()%3);
            const auto old=methods_mat(g,"D",1);
            const auto fresh=GWLD::Engine(g,GWLD::Method::D).matrix(GWLD::Method::D,1);
            std::vector<GWLD::Summary> summaries;
            for (arma::uword j=0;j<g.n_cols;++j) summaries.push_back(GWLD::summarize(g.colptr(j),g.n_rows));
            const std::string dataset=std::to_string(config[0])+"x"+std::to_string(config[1])+"_missing"+std::to_string(config[2]);
            for (int metric=0;metric<3;++metric) {
                double maximum=-1; std::size_t wi=0,wj=0;
                for (arma::uword j=1;j<g.n_cols;++j) for (arma::uword i=0;i<j;++i) {
                    const double p=summaries[i].p,q=summaries[j].p;
                    const auto value=[&](double d) {
                        if (metric==0) return d*d/(p*q*(1-p)*(1-q));
                        if (metric==1) return d;
                        return d/(d>0 ? std::min(p*(1-q),q*(1-p)) : std::max(-p*q,-(1-p)*(1-q)));
                    };
                    const double error=std::abs(value(old(i,j))-value(fresh(i,j)));
                    if (error>maximum) { maximum=error;wi=i;wj=j; }
                }
                emit(g,wi,wj,old(wi,wj),dataset,metric==0 ? "worst_r2" : metric==1 ? "worst_D" : "worst_Dprime");
            }
        }
        const auto table_case=[&](const std::array<int,9>& table,const std::string& name) {
            std::size_t n=0; for (int x:table) n+=x;
            arma::Mat<int> g(n,2);std::size_t k=0;
            for (int i=0;i<3;++i) for (int j=0;j<3;++j)
                for (int r=0;r<table[3*i+j];++r,++k) { g(k,0)=i;g(k,1)=j; }
            const auto old=methods_mat(g,"D",1);
            emit(g,0,1,old(0,1),name,"boundary_or_multimodal");
        };
        table_case({{0,2,0,2,262,5,2,5,2}},"multimodal_counterexample");
        table_case({{10,0,0,0,20,0,0,0,10}},"perfect_coupling");
        table_case({{0,0,10,0,20,0,10,0,0}},"perfect_repulsion");
        table_case({{0,0,0,0,20,0,0,0,0}},"tied_boundary_maxima");
        std::mt19937 random(654321);
        for (int trial=0;trial<24;++trial) {
            arma::Mat<int> g(256,2);
            for (auto& x:g) x=random()%8==0 ? -1 : static_cast<int>(random()%3);
            const auto old=methods_mat(g,"D",1);
            emit(g,0,1,old(0,1),"random_missing_"+std::to_string(trial+1),"seed654321");
        }
        cases.close();observations.close();
        std::cout<<"Exported "<<id<<" cases and raw genotypes\n";
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
