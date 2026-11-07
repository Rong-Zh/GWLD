#include <armadillo>
#include "gwld/engine.hpp"
#include <iostream>
#include <random>

void require(bool value,const char* message) { if (!value) throw std::runtime_error(message); }
void same(double a,double b) {
    require((std::isnan(a)&&std::isnan(b)) || (std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=1e-12),"Pair result differs from scalar reference");
}
void same_counts(const GWLD::Counts& a,const GWLD::Counts& b) {
    for (int i=0;i<9;++i) require(a.cells[i]==b.cells[i],"Contingency counts differ");
    for (int i=0;i<3;++i) require(a.rows[i]==b.rows[i] && a.cols[i]==b.cols[i],"Marginal counts differ");
    require(a.n==b.n,"Complete-pair count differs");
}
int main() {
    try {
        std::mt19937 random(2026);
        const int missing=std::numeric_limits<int>::min();
        for (std::size_t samples:{0,1,31,63,64,65,127,128,129,255,256,257,511,512,513,1023,1024,1025,4097,65535}) {
            arma::Mat<int> g(samples,7);
            for (auto& value:g) {
                const unsigned code=random()%5;
                value=code<3 ? static_cast<int>(code) : (code==3 ? -1 : missing);
            }
            GWLD::PackedGenotypes portable,native;
            portable.build(g,GWLD::count_packed_portable); native.build(g);
            std::vector<GWLD::Summary> summaries;
            for (arma::uword j=0;j<g.n_cols;++j) summaries.push_back(GWLD::summarize(samples ? g.colptr(j) : nullptr,samples));
            for (arma::uword j=1;j<g.n_cols;++j) for(arma::uword i=0;i<j;++i) {
                const auto scalar=GWLD::count(samples ? g.colptr(i) : nullptr,samples ? g.colptr(j) : nullptr,samples);
                if (samples) {
                    same_counts(portable.count(i,j),scalar); same_counts(native.count(i,j),scalar);
                }
            }
            for (auto method:{GWLD::Method::MI,GWLD::Method::RMI,GWLD::Method::R2,GWLD::Method::D,GWLD::Method::Dprime}) {
                const GWLD::Engine engine(g,method);
                const auto matrix=engine.matrix(method,1);
                for (arma::uword j=0;j<g.n_cols;++j) for (arma::uword i=0;i<g.n_cols;++i) {
                    if (i>=j) require(std::isnan(matrix(i,j)),"Lower triangle must stay missing");
                    else same(matrix(i,j),GWLD::pair(samples ? g.colptr(i) : nullptr,samples ? g.colptr(j) : nullptr,
                                                   samples,summaries[i],summaries[j],method));
                }
            }
        }
        arma::Mat<int> edge(513,16);
        for(auto& x:edge) x=static_cast<int>(random()%4)-1;
        edge.col(0).fill(-1); edge.col(1).fill(0);
        // A category observed only outside complete pairs must remain in RMI's summaries.
        edge.col(2).fill(0); edge.col(3).fill(1);
        edge(512,2)=2; edge(512,3)=missing;
        for (auto method:{GWLD::Method::MI,GWLD::Method::RMI,GWLD::Method::R2,GWLD::Method::D,GWLD::Method::Dprime}) {
            const GWLD::Engine engine(edge,method);
            const GWLD::Engine no_rmi_cache(edge,GWLD::Method::MI);
            const auto matrix=engine.matrix(method,1);
            #ifdef _OPENMP
            const auto parallel=engine.matrix(method,3);
            #endif
            for (arma::uword j=1;j<edge.n_cols;++j) for (arma::uword i=0;i<j;++i) {
                const auto sa=GWLD::summarize(edge.colptr(i),edge.n_rows),sb=GWLD::summarize(edge.colptr(j),edge.n_rows);
                same(matrix(i,j),GWLD::pair(edge.colptr(i),edge.colptr(j),edge.n_rows,sa,sb,method));
                same(matrix(i,j),no_rmi_cache.at(i,j,method));
                #ifdef _OPENMP
                same(matrix(i,j),parallel(i,j));
                #endif
            }
        }
        // Exact counts across NEON's 16-bit accumulator flush boundaries.
        for (std::size_t samples:{524159,524160,524161,1048321}) {
            arma::Mat<int> g(samples,4);
            for (auto& x:g) x=static_cast<int>(random()%3);
            const auto a=GWLD::summarize(g.colptr(0),samples),b=GWLD::summarize(g.colptr(1),samples);
            const auto reference=GWLD::count(g.colptr(0),g.colptr(1),samples);
            GWLD::PackedGenotypes portable,native;
            portable.build(g,GWLD::count_packed_portable); native.build(g);
            same_counts(reference,native.count(0,1));
            same_counts(reference,native.count_complete(0,1,a,b));
            same_counts(reference,portable.count_complete(0,1,a,b));
            g.fill(0); native.build(g);
            same_counts(GWLD::count(g.colptr(0),g.colptr(1),samples),native.count(0,1));
        }
        // Complete columns spanning every 1/2/3-category marginal combination.
        arma::Mat<int> categories(513,9);
        for (arma::uword j=0;j<categories.n_cols;++j) for (arma::uword i=0;i<categories.n_rows;++i)
            categories(i,j)=static_cast<int>((random()%(1+j%3)+j/3)%3);
        for (auto method:{GWLD::Method::MI,GWLD::Method::RMI}) {
            const GWLD::Engine engine(categories,method);
            const auto result=engine.matrix(method,1);
            for (arma::uword j=1;j<categories.n_cols;++j) for (arma::uword i=0;i<j;++i)
                same(result(i,j),GWLD::pair(categories.colptr(i),categories.colptr(j),categories.n_rows,
                     GWLD::summarize(categories.colptr(i),categories.n_rows),
                     GWLD::summarize(categories.colptr(j),categories.n_rows),method));
        }
        edge(0,0)=3;
        bool rejected=false;
        try { GWLD::PackedGenotypes invalid; invalid.build(edge); } catch(const std::invalid_argument&) {rejected=true;}
        require(rejected,"Invalid packed genotype accepted");
        std::cout<<"Portable/native packed counts, padding, missing values and scalar equivalence passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n'; return 1;}
}
