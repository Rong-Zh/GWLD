#include <iostream>
#include "ArgumentParser.h"
#include "dataset.hpp"
#include "gwld/engine.hpp"

namespace {
using namespace GWLD;
void matrix_file(const arma::Mat<double>& matrix,const std::string& path) {
    auto out=output(path);
    for (arma::uword i=0; i<matrix.n_rows; ++i) {
        for (arma::uword j=0; j<matrix.n_cols; ++j) {
            if (j) out<<'\t';
            out<<matrix(i,j);
        }
        out<<'\n';
    }
    out.close();
}
void pairs_file(const Dataset& data,const Engine& engine,Method method,int threads,
                bool decay,int maxdist,double threshold,const std::string& path,
                const arma::Mat<double>* matrix) {
    #ifndef _OPENMP
    (void)threads;
    #endif
    auto out=output(path);
    const auto n=data.sites.size();
    const Window window(data.sites);
    // Bounded batches keep pair output from allocating a V x V matrix.
    const std::size_t batch=64;
    for (std::size_t start=0; start<n; start+=batch) {
        const auto count=std::min(batch,n-start);
        std::vector<std::string> rows(count);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        #endif
        for (std::int64_t offset=0; offset<static_cast<std::int64_t>(count); ++offset) {
            const auto i=start+offset; const auto& a=data.sites[i];
            std::ostringstream buffer; buffer.precision(12);
            auto emit=[&](std::size_t j) {
                const auto& b=data.sites[j];
                if (!decay && a.chr==b.chr) return;
                const double value=matrix ? (*matrix)(i,j) : engine.at(i,j,method);
                if (decay) buffer<<a.id<<'\t'<<b.id<<'\t'<<std::llabs(b.pos-a.pos)<<'\t'<<value<<'\n';
                else if (value>threshold)
                    buffer<<a.chr<<'\t'<<a.pos<<'\t'<<a.id<<'\t'<<b.chr<<'\t'<<b.pos<<'\t'<<b.id<<'\t'<<value<<'\n';
            };
            if (decay) for (auto j:window.partners(i,a,static_cast<std::int64_t>(maxdist)*1000)) emit(j);
            else for (std::size_t j=i+1; j<n; ++j) emit(j);
            rows[offset]=buffer.str();
        }
        for (const auto& row:rows) out<<row;
    }
    out.close();
}
int run(int argc,char* argv[]) {
    ArgumentParser args("GWLD","Genome-wide linkage disequilibrium analysis");
    args.set_version("2.0.0"); args.add_version('V',"version","Show version");
    args.add_help('h',"help","Show help");
    args.add_optional<std::string>('b',"bfile","PLINK BED/BIM/FAM prefix","");
    args.add_optional<std::string>('f',"file","PLINK PED/MAP prefix","");
    args.add_optional<std::string>('v',"vcf","Biallelic diploid VCF, VCF.gz or BCF","");
    args.add_optional<std::string>('m',"method","RMI, MI, r2 (r^2), D, D'","RMI");
    args.add_optional<int>('p',"threads","Number of worker threads",1);
    args.add_optional<int>("max-dist","Maximum decay distance in kb",300);
    args.add_optional<double>("threshold","Circos threshold",0.1);
    args.add_optional<std::string>('o',"out-prefix","Output prefix","gwld");
    args.add_flag('d',"decay","Calculate within-chromosome decay pairs");
    args.add_flag('c',"circos","Calculate cross-chromosome pairs");
    args.add_flag("by-chrom","Write a heatmap for each chromosome");
    args.add_flag("heatmap","Also write the full heatmap with other outputs");
    args.add_flag("code012","Write numeric genotypes");
    args.add_flag("allele-freq","Write REF/ALT allele frequencies and missing rate");
    args.add_flag("2vcf","Write genotypes as VCF");
    if (argc==1) { args.print_help(); return 0; }
    if (!args.parse_args(argc,argv)) return 0;
    const auto bed=args.get<std::string>("bfile"),ped=args.get<std::string>("file"),vcf=args.get<std::string>("vcf");
    if ((!bed.empty())+(!ped.empty())+(!vcf.empty())!=1)
        throw std::invalid_argument("Specify exactly one input: --bfile, --file or --vcf");
    const int threads=args.get<int>("threads"),maxdist=args.get<int>("max-dist");
    if (threads<1 || maxdist<0) throw std::invalid_argument("Threads must be positive; max-dist must be nonnegative");
    const double threshold=args.get<double>("threshold");
    if (!std::isfinite(threshold)) throw std::invalid_argument("Threshold must be finite");
    #ifndef _OPENMP
    if (threads>1) throw std::invalid_argument("This build has no OpenMP support; use --threads 1");
    #endif
    const auto method=parse_method(args.get<std::string>("method"));
    const auto prefix=args.get<std::string>("out-prefix");
    auto data=!bed.empty() ? load_bed(bed) : (!ped.empty() ? load_ped(ped) : load_vcf(vcf,threads));
    Engine engine(data.genotypes);
    std::cerr<<"GWLD 2.0.0: "<<data.samples.size()<<" samples, "<<data.sites.size()<<" variants\n";
    const bool decay=args.get<bool>("decay"),circos=args.get<bool>("circos"),bychrom=args.get<bool>("by-chrom");
    const bool code=args.get<bool>("code012"),freq=args.get<bool>("allele-freq"),convert=args.get<bool>("2vcf");
    if (code || convert) write_genotypes(data,prefix+".vcf",code);
    if (freq) {
        auto out=output(prefix+".freq");
        for (arma::uword j=0; j<data.genotypes.n_cols; ++j) {
            const auto summary=summarize(data.genotypes.colptr(j),data.genotypes.n_rows);
            std::size_t missing=0;
            for (arma::uword i=0; i<data.genotypes.n_rows; ++i) missing+=data.genotypes(i,j)<0;
            out<<summary.p<<'\t'<<1-summary.p<<'\t'<<static_cast<double>(missing)/data.samples.size()<<'\n';
        }
        out.close();
    }
    arma::Mat<double> matrix;
    const bool full=args.get<bool>("heatmap") || !(decay || circos || bychrom || code || freq || convert);
    if (full) { matrix=engine.matrix(method,threads); matrix_file(matrix,prefix+".heatmap"); }
    if (decay) pairs_file(data,engine,method,threads,true,maxdist,threshold,prefix+".decay",full ? &matrix : nullptr);
    if (circos) pairs_file(data,engine,method,threads,false,maxdist,threshold,prefix+".circos",full ? &matrix : nullptr);
    if (bychrom) {
        std::map<std::string,std::vector<arma::uword>> chromosomes;
        for (std::size_t i=0; i<data.sites.size(); ++i) chromosomes[data.sites[i].chr].push_back(i);
        for (const auto& group:chromosomes) {
            if (group.first.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-")!=std::string::npos)
                throw std::invalid_argument("Chromosome name cannot be used in an output filename");
            arma::uvec index(group.second);
            arma::Mat<double> local;
            if (full) local=matrix.submat(index,index);
            else {
                arma::Mat<int> genotypes=data.genotypes.cols(index);
                local=Engine(genotypes).matrix(method,threads);
            }
            matrix_file(local,prefix+".chr"+group.first+".heatmap");
        }
    }
    return 0;
}
} // namespace
int main(int argc,char* argv[]) {
    try { return run(argc,argv); }
    catch (const std::exception& error) { std::cerr<<"GWLD: "<<error.what()<<'\n'; return 1; }
}
