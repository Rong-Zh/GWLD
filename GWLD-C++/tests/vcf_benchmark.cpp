#include "legacy/vcf_stream.hpp"
#include <chrono>
#include <iostream>

int main(int argc,char* argv[]) {
    if (argc!=2) { std::cerr<<"Usage: gwld_vcf_benchmark input.vcf\n"; return 1; }
    try {
        std::vector<double> old_times,new_times;
        for (int repeat=0;repeat<3;++repeat) {
            GWLD::Dataset old,fresh;
            auto measure=[&](bool legacy) {
                const auto start=std::chrono::steady_clock::now();
                if (legacy) old=GWLD::legacy::load_vcf(argv[1]);
                else fresh=GWLD::load_vcf(argv[1]);
                const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                (legacy ? old_times : new_times).push_back(elapsed);
            };
            measure(repeat%2==0); measure(repeat%2!=0);
            if (old.genotypes.n_rows!=fresh.genotypes.n_rows || old.genotypes.n_cols!=fresh.genotypes.n_cols ||
                old.samples!=fresh.samples || old.alleles!=fresh.alleles ||
                !std::equal(old.genotypes.begin(),old.genotypes.end(),fresh.genotypes.begin()))
                throw std::runtime_error("Genotype or metadata mismatch");
            for (std::size_t i=0;i<old.sites.size();++i)
                if (old.sites[i].chr!=fresh.sites[i].chr || old.sites[i].pos!=fresh.sites[i].pos || old.sites[i].id!=fresh.sites[i].id)
                    throw std::runtime_error("Site metadata mismatch");
        }
        std::sort(old_times.begin(),old_times.end()); std::sort(new_times.begin(),new_times.end());
        std::cout<<"median_stream_seconds\tmedian_htslib_seconds\tspeedup\n"
                 <<old_times[1]<<'\t'<<new_times[1]<<'\t'<<old_times[1]/new_times[1]<<'\n';
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
