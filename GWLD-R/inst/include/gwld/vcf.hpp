#pragma once
#include <armadillo>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <array>
#include <cstdlib>
#include <cstring>
#include <memory>
#include "core.hpp"

namespace GWLD {
struct VariantInfo {
    std::vector<Site> sites;
    std::vector<std::array<std::string,2>> alleles;
    std::vector<std::string> samples;
};
struct Dataset : VariantInfo { arma::Mat<int> genotypes; };
struct VcfData : VariantInfo { std::vector<std::int8_t> genotypes; };

// One reader for CLI and R: VCF, gzip/BGZF VCF and BCF. All HTSlib buffers
// are reused between records; the computation path never formats GT strings.
class VcfReader {
    struct Close { void operator()(htsFile* p) const { hts_close(p); } };
    struct Header { void operator()(bcf_hdr_t* p) const { bcf_hdr_destroy(p); } };
    struct Record { void operator()(bcf1_t* p) const { bcf_destroy(p); } };
    struct Buffers {
        int32_t* gt=nullptr;
        int capacity=0;
        kstring_t text={0,0,nullptr};
        kstring_t line={0,0,nullptr};
        ~Buffers() { std::free(gt); std::free(text.s); std::free(line.s); }
    } buffers;
    std::unique_ptr<htsFile,Close> file;
    std::unique_ptr<bcf_hdr_t,Header> header;
    std::unique_ptr<bcf1_t,Record> record;
    std::vector<std::int8_t> genotype_codes;
    std::size_t ordinal=0;
    int stride=0;
    static bool missing(int32_t value) { return value==bcf_int32_missing || bcf_gt_is_missing(value); }

    void fail(const char* message) const {
        throw std::runtime_error("VCF/BCF record "+std::to_string(ordinal)+": "+message);
    }
    void decode() {
        if (record->n_allele>2) fail("only biallelic variants are supported");
        const int count=bcf_get_genotypes(header.get(),record.get(),&buffers.gt,&buffers.capacity);
        const auto n=sample_count();
        if (count<=0 || static_cast<std::size_t>(count)%n) fail("missing or invalid FORMAT/GT");
        stride=static_cast<int>(count/n);
        for (std::size_t i=0;i<n;++i) {
            const int32_t* gt=buffers.gt+i*stride;
            if (stride>2 && gt[2]!=bcf_int32_vector_end) fail("only diploid genotypes are supported");
            if (stride<2 || gt[1]==bcf_int32_vector_end) {
                if (gt[0]==bcf_int32_vector_end || !missing(gt[0])) fail("haploid GT is not supported");
                genotype_codes[i]=-1; continue; // a single '.' is an entirely missing sample
            }
            if (gt[0]==bcf_int32_vector_end) fail("invalid GT ploidy");
            const int a=bcf_gt_allele(gt[0]),b=bcf_gt_allele(gt[1]);
            if ((!missing(gt[0]) && (a<0 || a>=record->n_allele)) ||
                (!missing(gt[1]) && (b<0 || b>=record->n_allele))) fail("GT allele index exceeds REF/ALT");
            genotype_codes[i]=(missing(gt[0]) || missing(gt[1])) ? -1 : a+b;
        }
    }
public:
    explicit VcfReader(const std::string& path,int threads=1) : file(hts_open(path.c_str(),"r")) {
        if (!file) throw std::runtime_error("Cannot open VCF/BCF: "+path);
        if (threads<1) throw std::invalid_argument("Threads must be positive");
        // Large reads matter especially on network filesystems and WSL mounts.
        if (hts_set_opt(file.get(),HTS_OPT_BLOCK_SIZE,4*1024*1024)<0)
            throw std::runtime_error("Cannot configure HTSlib input buffer");
        // The calling thread plus N-1 BGZF workers stay within the requested budget.
        if (threads>1 && hts_set_threads(file.get(),threads-1)<0)
            throw std::runtime_error("Cannot enable HTSlib decompression threads");
        header.reset(bcf_hdr_read(file.get()));
        if (!header || bcf_hdr_nsamples(header.get())<=0)
            throw std::runtime_error("VCF/BCF has no valid sample header");
        // Support older VCFs that omitted the standard GT declaration.
        const int gt=bcf_hdr_id2int(header.get(),BCF_DT_ID,"GT");
        if (!bcf_hdr_idinfo_exists(header.get(),BCF_HL_FMT,gt)) {
            if (bcf_hdr_append(header.get(),"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")<0 ||
                bcf_hdr_sync(header.get())<0) throw std::runtime_error("Cannot initialize GT header");
        }
        record.reset(bcf_init());
        if (!record) throw std::bad_alloc();
        genotype_codes.resize(sample_count());
    }
    std::size_t sample_count() const { return static_cast<std::size_t>(bcf_hdr_nsamples(header.get())); }
    std::vector<std::string> sample_names() const {
        return std::vector<std::string>(header->samples,header->samples+sample_count());
    }
    bool next(bool genotypes=true) {
        int status;
        if (hts_get_format(file.get())->format==vcf) {
            // Check columns before HTSlib tokenizes the reusable line buffer.
            // vcf_parse can otherwise silently discard surplus samples.
            do { status=hts_getline(file.get(),'\n',&buffers.line); }
            while (status>=0 && !buffers.line.l);
            if (status==-1) return false;
            ++ordinal;
            if (status<0) fail("malformed or truncated input");
            std::size_t tabs=0;
            for (std::size_t i=0;i<buffers.line.l;++i) tabs+=buffers.line.s[i]=='\t';
            if (tabs!=sample_count()+8) fail("sample count differs from header");
            bcf_clear(record.get());
            status=vcf_parse(&buffers.line,header.get(),record.get());
        } else {
            status=bcf_read(file.get(),header.get(),record.get());
            if (status==-1) return false;
            ++ordinal;
        }
        if (status<0 || (record->errcode & ~(BCF_ERR_CTG_UNDEF|BCF_ERR_TAG_UNDEF))) fail("malformed or truncated input");
        if (record->n_sample!=sample_count()) fail("sample count differs from header");
        if (bcf_unpack(record.get(),BCF_UN_STR)<0) fail("cannot decode variant metadata");
        if (record->pos<0 || record->pos>=std::numeric_limits<std::int64_t>::max()/2) fail("invalid site position");
        if (!bcf_hdr_id2name(header.get(),record->rid)) fail("invalid chromosome");
        if (genotypes) decode();
        return true;
    }
    Site site() const { return {bcf_hdr_id2name(header.get(),record->rid),record->d.id,record->pos+1}; }
    std::array<std::string,2> alleles() const {
        return {{record->d.allele[0],record->n_allele>1 ? record->d.allele[1] : "."}};
    }
    const std::vector<std::int8_t>& genotypes() const { return genotype_codes; }
    int allele(std::size_t sample,int copy) const {
        if (copy>=stride) return -1;
        const auto value=buffers.gt[sample*stride+copy];
        return value==bcf_int32_vector_end || missing(value) ? -1 : bcf_gt_allele(value);
    }
    // Only the legacy R character-table adapter requests formatted fields.
    std::vector<std::string> text_fields() {
        buffers.text.l=0;
        if (vcf_format(header.get(),record.get(),&buffers.text)<0) fail("cannot format record");
        const char* begin=buffers.text.s; const char* end=begin+buffers.text.l;
        while (end>begin && (end[-1]=='\n' || end[-1]=='\r')) --end;
        std::vector<std::string> fields; fields.reserve(9+sample_count());
        for (;;) {
            const auto tab=static_cast<const char*>(std::memchr(begin,'\t',end-begin));
            if (!tab) { fields.emplace_back(begin,end); break; }
            fields.emplace_back(begin,tab); begin=tab+1;
        }
        return fields;
    }
};

inline VcfData read_vcf_data(const std::string& path,int threads=1) {
    VcfReader reader(path,threads);
    VcfData data; data.samples=reader.sample_names();
    // One byte per pending genotype limits staging memory. No two-pass read,
    // per-sample allocations, or four-byte matrix reallocations while streaming.
    while (reader.next()) {
        data.sites.push_back(reader.site()); data.alleles.push_back(reader.alleles());
        const auto& genotypes=reader.genotypes();
        data.genotypes.insert(data.genotypes.end(),genotypes.begin(),genotypes.end());
    }
    if (data.sites.empty()) throw std::runtime_error("VCF/BCF has no variants");
    const auto samples=data.samples.size(),variants=data.sites.size();
    if (variants>std::numeric_limits<arma::uword>::max()/samples ||
        variants>std::numeric_limits<std::size_t>::max()/sizeof(int)/samples)
        throw std::length_error("Genotype matrix is too large");
    return data;
}
inline Dataset load_vcf(const std::string& path,int threads=1) {
    auto compact=read_vcf_data(path,threads);
    Dataset data;
    data.sites=std::move(compact.sites); data.alleles=std::move(compact.alleles); data.samples=std::move(compact.samples);
    data.genotypes.set_size(data.samples.size(),data.sites.size());
    std::copy(compact.genotypes.begin(),compact.genotypes.end(),data.genotypes.memptr());
    return data;
}
} // namespace GWLD
