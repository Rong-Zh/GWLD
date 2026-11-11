#include <armadillo>
#include "gwld/engine.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <numeric>

// Deterministic, end-to-end engine timing, including summaries, caches and packing.
int main(int argc,char** argv) {
    try {
        if (argc<5) throw std::invalid_argument("Usage: kernel_benchmark samples variants missing_percent threads [output_prefix] [uniform|hwe] [constant_percent]");
        const std::size_t samples=std::stoul(argv[1]),variants=std::stoul(argv[2]);
        const int missing=std::stoi(argv[3]),threads=std::stoi(argv[4]);
        if (missing<0 || missing>100 || threads<1) throw std::invalid_argument("Invalid missing percentage or threads");
        const std::string distribution=argc>6 ? argv[6] : "uniform";
        if (distribution!="uniform" && distribution!="hwe") throw std::invalid_argument("Choose uniform or hwe genotypes");
        const int constant=argc>7 ? std::stoi(argv[7]) : 0;
        if (constant<0 || constant>100) throw std::invalid_argument("Constant percentage must be between 0 and 100");
        arma::Mat<int> g(samples,variants); std::mt19937 random(2026);
        for (std::size_t j=0;j<variants;++j) {
            // Different minor-allele frequencies across variants, under HWE.
            const double p=0.01+0.49*(variants>1 ? static_cast<double>(j)/(variants-1) : 0.5);
            const double cutoff=p*static_cast<double>(std::mt19937::max());
            for (std::size_t i=0;i<samples;++i) {
                if (random()%100<static_cast<unsigned>(missing)) g(i,j)=-1;
                else if (distribution=="uniform") g(i,j)=static_cast<int>(random()%3);
                else g(i,j)=static_cast<int>(random()<cutoff)+static_cast<int>(random()<cutoff);
            }
        }
        // Constant columns are scattered with a separate seed so input order
        // does not artificially favor one branch direction. Category 1 is an
        // all-heterozygous column, not an allele-fixed locus.
        if (constant) {
            std::vector<std::size_t> order(variants);
            std::iota(order.begin(),order.end(),0);
            std::mt19937 positions(2027);
            std::shuffle(order.begin(),order.end(),positions);
            const auto count=(variants/100)*constant+(variants%100)*constant/100;
            for (std::size_t j=0;j<count;++j) g.col(order[j]).fill(static_cast<int>(j%3));
        }
        std::cout<<"method\tsamples\tvariants\tmissing_percent\tconstant_percent\tthreads\tmedian_seconds\tchecksum\n"<<std::setprecision(17);
        for(const std::string name:{"MI","RMI","r2","D","D'"}) {
            const auto method=GWLD::parse_method(name);
            std::vector<double> times; arma::Mat<double> result;
            for(int repeat=0;repeat<4;++repeat) {
                const auto start=std::chrono::steady_clock::now();
                result=GWLD::Engine(g,method).matrix(method,threads);
                const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                if(repeat) times.push_back(elapsed); // warm caches and OpenMP runtime
            }
            std::sort(times.begin(),times.end());
            double checksum=0;
            for(arma::uword j=1;j<result.n_cols;++j) for(arma::uword i=0;i<j;++i)
                if(std::isfinite(result(i,j))) checksum+=result(i,j);
            std::cout<<name<<'\t'<<samples<<'\t'<<variants<<'\t'<<missing<<'\t'<<constant<<'\t'<<threads<<'\t'<<times[1]<<'\t'<<checksum<<std::endl;
            if(argc>5 && std::string(argv[5])!="-") {
                std::ofstream out(std::string(argv[5])+"-"+name+".bin",std::ios::binary);
                out.write(reinterpret_cast<const char*>(result.memptr()),result.n_elem*sizeof(double));
                if(!out) throw std::runtime_error("Cannot write benchmark matrix");
            }
        }
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
