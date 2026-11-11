#include "gwld/packed.hpp"
#include <iostream>
#include <random>
#include <utility>

struct Matrix {
    std::size_t n_rows,n_cols=4;
    std::vector<int> values;
    explicit Matrix(std::size_t rows) : n_rows(rows),values(rows*n_cols) {}
    const int* colptr(std::size_t column) const { return values.data()+column*n_rows; }
};
void equal_counts(const GWLD::Counts& a,const GWLD::Counts& b) {
    for (unsigned k=0;k<9;++k) if (a.cells[k]!=b.cells[k]) throw std::runtime_error("SIMD cell count mismatch");
    for (unsigned k=0;k<3;++k) if (a.rows[k]!=b.rows[k] || a.cols[k]!=b.cols[k]) throw std::runtime_error("SIMD marginal mismatch");
    if (a.n!=b.n) throw std::runtime_error("SIMD sample count mismatch");
}
int main() {
    try {
        std::vector<std::pair<const char*,GWLD::CountKernel>> kernels={{"portable",GWLD::count_packed_portable},{"automatic",nullptr}};
        #if defined(GWLD_X86_SIMD)
        // A hypervisor can mask POPCNT independently from the vector features.
        // Check pointer selection without executing an unsupported kernel.
        for (std::size_t words:{0,8,16,128}) {
            if (GWLD::select_x86_count_kernel(words,false,true,true)!=GWLD::count_packed_portable)
                throw std::runtime_error("SIMD dispatch requires scalar POPCNT");
            if (GWLD::select_x86_count_kernel(words,true,false,true)!=GWLD::count_packed_native)
                throw std::runtime_error("Scalar POPCNT dispatch mismatch");
        }
        if (GWLD::select_x86_count_kernel(15,true,true,false)!=GWLD::count_packed_native ||
            GWLD::select_x86_count_kernel(16,true,true,false)!=GWLD::count_packed_avx2<>)
            throw std::runtime_error("AVX2 dispatch size boundary mismatch");
        #if defined(GWLD_X86_AVX512)
        if (GWLD::select_x86_count_kernel(7,true,true,true)!=GWLD::count_packed_native ||
            GWLD::select_x86_count_kernel(8,true,true,true)!=GWLD::count_packed_avx512<>)
            throw std::runtime_error("AVX512 dispatch size boundary mismatch");
        #endif
        if (__builtin_cpu_supports("popcnt")) kernels.emplace_back("popcnt",GWLD::count_packed_native);
        if (__builtin_cpu_supports("popcnt") && __builtin_cpu_supports("avx2")) kernels.emplace_back("avx2",GWLD::count_packed_avx2<>);
        #if defined(GWLD_X86_AVX512)
        if (__builtin_cpu_supports("popcnt") && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vpopcntdq"))
            kernels.emplace_back("avx512vpopcntdq",GWLD::count_packed_avx512<>);
        #endif
        #elif defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
        kernels.emplace_back("native",GWLD::count_packed_native);
        kernels.emplace_back("neon",GWLD::count_packed_neon<>);
        #endif
        std::mt19937 random(2026);
        for (std::size_t n:{0,1,63,64,65,127,128,129,255,256,257,511,512,513,1023,1024,1025,
                           7935,7936,7937,8191,8192,8193,524159,524160,524161}) {
            for (int pattern=0;pattern<3;++pattern) {
                Matrix g(n);
                for (auto& x:g.values) x=pattern==2 ? 0 : static_cast<int>(random()%(pattern==1 ? 4:3));
                if (pattern==1) for (auto& x:g.values) if (x==3) x=-1;
                std::vector<GWLD::Summary> summaries;
                for (unsigned j=0;j<4;++j) summaries.push_back(GWLD::summarize(n ? g.colptr(j):nullptr,n));
                for (const auto& kernel:kernels) {
                    GWLD::PackedGenotypes packed; packed.build(g,kernel.second);
                    if (!n) continue;
                    for (unsigned j=1;j<4;++j) for (unsigned i=0;i<j;++i) {
                        const auto reference=GWLD::count(g.colptr(i),g.colptr(j),n);
                        equal_counts(reference,packed.count(i,j));
                        if (pattern!=1) equal_counts(reference,packed.count_complete(i,j,summaries[i],summaries[j]));
                    }
                }
            }
        }
        // Inspect every bit, not just pair tables: a shared sample permutation
        // could otherwise pass contingency-count comparisons unnoticed.
        for (std::size_t n:{1,15,16,17,31,32,33,63,64,65,127,128,129}) {
            Matrix g(n);
            for (auto& x:g.values) {
                const auto code=random()%5;
                x=code<3 ? static_cast<int>(code) : code==3 ? -1 : std::numeric_limits<int>::min();
            }
            const auto words=(n+63)/64;
            std::vector<std::uint64_t> bits(words*3);
            GWLD::packing_detail::column(g.colptr(0),n,words,bits.data());
            for (std::size_t k=0;k<words*64;++k) for (unsigned code=0;code<3;++code) {
                const bool actual=(bits[code*words+k/64]>>(k%64))&1;
                const bool expected=k<n && g.values[k]==static_cast<int>(code);
                if (actual!=expected) throw std::runtime_error("Packed bit position or tail mismatch");
            }
            // Invalid codes must be rejected in every vector lane and tail.
            for (int code:{3,-2,std::numeric_limits<int>::min()+1,std::numeric_limits<int>::max()})
                for (std::size_t k=0;k<n;++k) {
                    const auto original=g.values[k]; g.values[k]=code;
                    bool rejected=false;
                    try { GWLD::packing_detail::column(g.colptr(0),n,words,bits.data()); }
                    catch (const std::invalid_argument&) { rejected=true; }
                    g.values[k]=original;
                    if (!rejected) throw std::runtime_error("Packing accepted an invalid genotype");
                }
        }
        // Reused storage must overwrite old planes, including a shorter tail.
        Matrix reused(129);
        GWLD::PackedGenotypes packed;
        packed.build(reused);
        std::fill(reused.values.begin(),reused.values.end(),-1);
        packed.build(reused);
        equal_counts(GWLD::Counts{},packed.count(0,1));
        reused=Matrix(65);
        for (auto& x:reused.values) x=static_cast<int>(random()%3);
        packed.build(reused);
        equal_counts(GWLD::count(reused.colptr(0),reused.colptr(1),65),packed.count(0,1));
        std::cout<<"Validated kernels:";
        for (const auto& kernel:kernels) std::cout<<' '<<kernel.first;
        std::cout<<"; exact counts, missing values, tails and accumulator boundaries passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
