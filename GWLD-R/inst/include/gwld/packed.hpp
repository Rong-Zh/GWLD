#pragma once
#include "core.hpp"
#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define GWLD_FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define GWLD_FORCE_INLINE __forceinline
#else
#define GWLD_FORCE_INLINE inline
#endif

namespace GWLD {
template<bool Native>
GWLD_FORCE_INLINE unsigned population_count(std::uint64_t value) {
    #if defined(__GNUC__) || defined(__clang__)
    if (Native) return static_cast<unsigned>(__builtin_popcountll(value));
    #elif defined(_MSC_VER) && defined(_M_X64)
    if (Native) return static_cast<unsigned>(__popcnt64(value));
    #endif
    // Portable SWAR fallback, including machines without POPCNT.
    value-=(value>>1)&UINT64_C(0x5555555555555555);
    value=(value&UINT64_C(0x3333333333333333))+((value>>2)&UINT64_C(0x3333333333333333));
    value=(value+(value>>4))&UINT64_C(0x0f0f0f0f0f0f0f0f);
    return static_cast<unsigned>((value*UINT64_C(0x0101010101010101))>>56);
}

template<bool Native>
GWLD_FORCE_INLINE Counts packed_count_kernel(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    std::uint64_t totals[9]={};
    for (std::size_t k=0;k<words;++k) {
        const auto a0=a[k],a1=a[words+k],a2=a[2*words+k];
        const auto b0=b[k],b1=b[words+k],b2=b[2*words+k];
        totals[0]+=population_count<Native>(a0&b0);
        totals[1]+=population_count<Native>(a0&b1);
        totals[2]+=population_count<Native>(a0&b2);
        totals[3]+=population_count<Native>(a1&b0);
        totals[4]+=population_count<Native>(a1&b1);
        totals[5]+=population_count<Native>(a1&b2);
        totals[6]+=population_count<Native>(a2&b0);
        totals[7]+=population_count<Native>(a2&b1);
        totals[8]+=population_count<Native>(a2&b2);
    }
    Counts result;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        const double value=static_cast<double>(totals[3*i+j]);
        result.cells[3*i+j]=value; result.rows[i]+=value; result.cols[j]+=value; result.n+=value;
    }
    return result;
}
inline Counts count_packed_portable(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    return packed_count_kernel<false>(a,b,words);
}
#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
__attribute__((target("popcnt")))
#endif
inline Counts count_packed_native(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    return packed_count_kernel<true>(a,b,words);
}
typedef Counts (*CountKernel)(const std::uint64_t*,const std::uint64_t*,std::size_t);
inline CountKernel select_count_kernel() {
    #if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
    return __builtin_cpu_supports("popcnt") ? count_packed_native : count_packed_portable;
    #elif defined(_MSC_VER) && defined(_M_X64)
    int registers[4]; __cpuid(registers,1);
    return (registers[2]&(1<<23)) ? count_packed_native : count_packed_portable;
    #elif defined(__GNUC__) || defined(__clang__)
    return count_packed_native; // compiler lowers the builtin for the target ISA
    #else
    return count_packed_portable;
    #endif
}

class PackedGenotypes {
    std::size_t words=0;
    std::vector<std::uint64_t> bits;
    CountKernel kernel=count_packed_portable;
public:
    template<class Matrix>
    void build(const Matrix& matrix,CountKernel selected=select_count_kernel()) {
        words=matrix.n_rows/64+(matrix.n_rows%64!=0);
        if (!words || !matrix.n_cols) { bits.clear(); return; }
        if (words>bits.max_size()/3 || matrix.n_cols>bits.max_size()/(3*words))
            throw std::length_error("Packed genotype matrix is too large");
        bits.assign(matrix.n_cols*3*words,0); kernel=selected;
        for (std::size_t j=0;j<matrix.n_cols;++j) {
            const int* column=matrix.colptr(j);
            for (std::size_t i=0;i<matrix.n_rows;++i) {
                const int code=column[i];
                if (code>=0 && code<=2) bits[j*3*words+code*words+i/64]|=UINT64_C(1)<<(i%64);
                else if (code!=-1 && code!=std::numeric_limits<int>::min())
                    throw std::invalid_argument("Genotypes must be 0, 1, 2 or missing (-1/NA)");
            }
        }
    }
    bool ready() const { return !bits.empty(); }
    Counts count(std::size_t i,std::size_t j) const {
        return kernel(bits.data()+i*3*words,bits.data()+j*3*words,words);
    }
};
} // namespace GWLD
#undef GWLD_FORCE_INLINE
