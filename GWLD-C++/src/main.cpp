#include <iostream>
#include "SubCommandParser.h"
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
    #if !defined(_OPENMP) || defined(GWLD_NO_OPENMP)
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
        #if defined(_OPENMP) && !defined(GWLD_NO_OPENMP)
        #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        #endif
        for (std::int64_t offset=0; offset<static_cast<std::int64_t>(count); ++offset) {
            const auto i=start+offset; const auto& a=data.sites[i];
            std::ostringstream buffer; buffer.precision(std::numeric_limits<double>::max_digits10);
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
enum class Command { Compute, Convert, Frequency };
struct Options {
    Command command;
    std::string bed,ped,vcf,prefix,format="vcf";
    Method method=Method::RMI;
    int threads=1,max_distance_kb=300;
    double circos_threshold=0.1;
    bool decay=false,circos=false,by_chromosome=false,heatmap=false;
};
Method cli_method(const std::string& name) {
    if (name=="rmi") return Method::RMI;
    if (name=="mi") return Method::MI;
    if (name=="r2") return Method::R2;
    if (name=="d") return Method::D;
    if (name=="dprime") return Method::Dprime;
    throw std::invalid_argument("Unknown method: "+name+"; choose rmi, mi, r2, d or dprime");
}
void common_arguments(ArgumentOptions& args) {
    args.add_help('h',"help","Show help");
    args.add_optional<std::string>('b',"bed-prefix","PLINK BED/BIM/FAM prefix (choose one input)","");
    args.add_optional<std::string>('p',"ped-prefix","PLINK PED/MAP prefix (choose one input)","");
    args.add_optional<std::string>('v',"vcf","Biallelic diploid VCF, VCF.gz or BCF (choose one input)","");
    args.add_optional<int>('t',"threads","Worker/decompression threads",1);
    args.add_optional<std::string>('o',"output-prefix","Output filename prefix","gwld");
}
Options common_options(ArgumentOptions& args,Command command) {
    if (args.was_set("bed-prefix")+args.was_set("ped-prefix")+args.was_set("vcf")!=1)
        throw std::invalid_argument("Specify exactly one input: --bed-prefix, --ped-prefix or --vcf");
    Options options;
    options.command=command;
    options.bed=args.get<std::string>("bed-prefix");
    options.ped=args.get<std::string>("ped-prefix");
    options.vcf=args.get<std::string>("vcf");
    options.threads=args.get<int>("threads");
    options.prefix=args.get<std::string>("output-prefix");
    return options;
}
void frequency_file(const Dataset& data,const std::string& path) {
    auto out=output(path);
    out<<"CHROM\tPOS\tID\tREF\tALT\tREF_FREQ\tALT_FREQ\tMISSING_RATE\n";
    for (arma::uword j=0; j<data.genotypes.n_cols; ++j) {
        const int* column=data.genotypes.colptr(j);
        const auto summary=GWLD::summarize(column,data.genotypes.n_rows);
        const auto called=static_cast<std::size_t>(summary.total);
        const double p=summary.p;
        const auto& site=data.sites[j];
        out<<site.chr<<'\t'<<site.pos<<'\t'<<site.id<<'\t'<<data.alleles[j][0]<<'\t'<<data.alleles[j][1]
           <<'\t'<<p<<'\t'<<1-p<<'\t'<<static_cast<double>(data.samples.size()-called)/data.samples.size()<<'\n';
    }
    out.close();
}
int execute(const Options& options) {
    const auto& bed=options.bed;
    const auto& ped=options.ped;
    const auto& vcf=options.vcf;
    if ((!bed.empty())+(!ped.empty())+(!vcf.empty())!=1)
        throw std::invalid_argument("Specify exactly one input: --bed-prefix, --ped-prefix or --vcf");
    const int threads=options.threads,maxdist=options.max_distance_kb;
    if (threads<1) throw std::invalid_argument("--threads must be positive");
    if (maxdist<0) throw std::invalid_argument("--max-distance-kb must be nonnegative");
    const double threshold=options.circos_threshold;
    if (!std::isfinite(threshold)) throw std::invalid_argument("--circos-threshold must be finite");
    if (options.prefix.empty()) throw std::invalid_argument("--output-prefix must not be empty");
    if (options.command==Command::Convert && options.format!="vcf" && options.format!="genotype")
        throw std::invalid_argument("--format must be vcf or genotype");
    #if !defined(_OPENMP) || defined(GWLD_NO_OPENMP)
    if (threads>1) throw std::invalid_argument("This build has no OpenMP support; use --threads 1");
    #endif
    const auto method=options.method;
    const auto& prefix=options.prefix;
    auto data=!bed.empty() ? load_bed(bed) : (!ped.empty() ? load_ped(ped) : load_vcf(vcf,threads));
    std::cerr<<"GWLD 2.0.0: "<<data.samples.size()<<" samples, "<<data.sites.size()<<" variants\n";
    if (options.command==Command::Convert) {
        if (options.format=="vcf") write_vcf(data,prefix+".vcf");
        else write_genotype(data,prefix+".genotype.tsv");
        return 0;
    }
    if (options.command==Command::Frequency) { frequency_file(data,prefix+".freq.tsv"); return 0; }
    const bool decay=options.decay,circos=options.circos,bychrom=options.by_chromosome;
    arma::Mat<double> matrix;
    const bool full=options.heatmap || !(decay || circos || bychrom);
    if (full || decay || circos) {
        const Engine engine(data.genotypes,method);
        if (full) { matrix=engine.matrix(method,threads); matrix_file(matrix,prefix+".heatmap"); }
        if (decay) pairs_file(data,engine,method,threads,true,maxdist,threshold,prefix+".decay",full ? &matrix : nullptr);
        if (circos) pairs_file(data,engine,method,threads,false,maxdist,threshold,prefix+".circos",full ? &matrix : nullptr);
    }
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
                local=Engine(genotypes,method).matrix(method,threads);
            }
            matrix_file(local,prefix+".chr"+group.first+".heatmap");
        }
    }
    return 0;
}
int run(int argc,char* argv[]) {
    SubCommandParser parser("GWLD","Genome-wide linkage disequilibrium analysis");
    parser.set_version("2.0.0");
    parser.add_version('V',"version","Show version");
    auto& compute=parser.add_subcommand("calculate","Calculate LD matrices and pairs",[](ArgumentOptions& args) {
        auto options=common_options(args,Command::Compute);
        options.method=cli_method(args.get<std::string>("method"));
        options.decay=args.get<bool>("decay");
        options.circos=args.get<bool>("circos");
        options.by_chromosome=args.get<bool>("by-chromosome");
        options.heatmap=args.get<bool>("heatmap");
        options.max_distance_kb=args.get<int>("max-distance-kb");
        options.circos_threshold=args.get<double>("circos-threshold");
        if (args.was_set("max-distance-kb") && !options.decay)
            throw std::invalid_argument("--max-distance-kb requires --decay");
        if (args.was_set("circos-threshold") && !options.circos)
            throw std::invalid_argument("--circos-threshold requires --circos");
        return execute(options);
    });
    common_arguments(compute);
    compute.add_optional<std::string>('m',"method","LD measure: rmi, mi, r2, d, dprime","rmi");
    compute.add_flag('d',"decay","Write within-chromosome pairs");
    compute.add_optional<int>("max-distance-kb","Maximum decay distance in kb (requires --decay)",300);
    compute.add_flag('c',"circos","Write cross-chromosome pairs");
    compute.add_optional<double>("circos-threshold","Minimum exclusive circos value (requires --circos)",0.1);
    compute.add_flag("by-chromosome","Write a heatmap for each chromosome");
    compute.add_flag("heatmap","Write the full LD matrix (default when no other LD output is selected)");
    auto& convert=parser.add_subcommand("convert","Convert genotypes to VCF or genotype TSV",[](ArgumentOptions& args) {
        auto options=common_options(args,Command::Convert);
        options.format=args.get<std::string>("format");
        return execute(options);
    });
    common_arguments(convert);
    convert.add_optional<std::string>('f',"format","Output format: vcf or genotype","vcf");
    auto& frequency=parser.add_subcommand("summarize","Write allele frequencies and missing rates",[](ArgumentOptions& args) {
        return execute(common_options(args,Command::Frequency));
    });
    common_arguments(frequency);
    parser.footnote("Examples:\n  GWLD calculate --bed-prefix example/example --method rmi --threads 4 -o result\n"
                    "  GWLD convert --vcf input.vcf.gz --format genotype -o result\n"
                    "  GWLD summarize --ped-prefix example/example-2 -o result\n"
                    "Use GWLD <command> --help for command options.");
    if (!parser.parse_args(argc,argv)) return 0;
    return parser.execute();
}
} // namespace
int main(int argc,char* argv[]) {
    try { return run(argc,argv); }
    catch (const std::exception& error) { std::cerr<<"GWLD: "<<error.what()<<'\n'; return 1; }
}
