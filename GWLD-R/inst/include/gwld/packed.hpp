#pragma once
#include "core.hpp"
#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#define GWLD_X86_SIMD 1
#include <immintrin.h>
#if (defined(__clang__) && __clang_major__>=7) || (!defined(__clang__) && __GNUC__>=8)
#define GWLD_X86_AVX512 1
#endif
#endif
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
#include <arm_neon.h>
#endif
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

// Complete kernels produce cells only; count_complete restores known marginals.
template<bool Native,bool Complete=false>
GWLD_FORCE_INLINE Counts packed_count_kernel(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    std::uint64_t totals[9]={};
    for (std::size_t k=0;k<words;++k) {
        const auto a0=a[k],a1=a[words+k],a2=a[2*words+k];
        const auto b0=b[k],b1=b[words+k],b2=b[2*words+k];
        totals[0]+=population_count<Native>(a0&b0);
        totals[1]+=population_count<Native>(a0&b1);
        if (!Complete) totals[2]+=population_count<Native>(a0&b2);
        totals[3]+=population_count<Native>(a1&b0);
        totals[4]+=population_count<Native>(a1&b1);
        if (!Complete) totals[5]+=population_count<Native>(a1&b2);
        if (!Complete) totals[6]+=population_count<Native>(a2&b0);
        if (!Complete) totals[7]+=population_count<Native>(a2&b1);
        if (!Complete) totals[8]+=population_count<Native>(a2&b2);
    }
    Counts result;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        const double value=static_cast<double>(totals[3*i+j]);
        result.cells[3*i+j]=value;
        if (!Complete) { result.rows[i]+=value; result.cols[j]+=value; result.n+=value; }
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
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
template<bool Complete=false>
GWLD_FORCE_INLINE Counts count_packed_neon(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    std::uint64_t totals[9]={};
    std::size_t k=0;
    while (k+1<words) {
        // Each lane adds at most 16 per pair of words. Bound chunks so that
        // uint16 accumulators cannot overflow even on cohorts above 500k samples.
        const auto end=k+2*std::min<std::size_t>((words-k)/2,4095);
        uint16x8_t partial[9]={};
        for (;k<end;k+=2) {
            const auto a0=vreinterpretq_u8_u64(vld1q_u64(a+k));
            const auto a1=vreinterpretq_u8_u64(vld1q_u64(a+words+k));
            const auto a2=vreinterpretq_u8_u64(vld1q_u64(a+2*words+k));
            const auto b0=vreinterpretq_u8_u64(vld1q_u64(b+k));
            const auto b1=vreinterpretq_u8_u64(vld1q_u64(b+words+k));
            const auto b2=vreinterpretq_u8_u64(vld1q_u64(b+2*words+k));
            partial[0]=vpadalq_u8(partial[0],vcntq_u8(vandq_u8(a0,b0)));
            partial[1]=vpadalq_u8(partial[1],vcntq_u8(vandq_u8(a0,b1)));
            if (!Complete) partial[2]=vpadalq_u8(partial[2],vcntq_u8(vandq_u8(a0,b2)));
            partial[3]=vpadalq_u8(partial[3],vcntq_u8(vandq_u8(a1,b0)));
            partial[4]=vpadalq_u8(partial[4],vcntq_u8(vandq_u8(a1,b1)));
            if (!Complete) partial[5]=vpadalq_u8(partial[5],vcntq_u8(vandq_u8(a1,b2)));
            if (!Complete) partial[6]=vpadalq_u8(partial[6],vcntq_u8(vandq_u8(a2,b0)));
            if (!Complete) partial[7]=vpadalq_u8(partial[7],vcntq_u8(vandq_u8(a2,b1)));
            if (!Complete) partial[8]=vpadalq_u8(partial[8],vcntq_u8(vandq_u8(a2,b2)));
        }
        for (unsigned cell=0;cell<9;++cell) totals[cell]+=vaddlvq_u16(partial[cell]);
    }
    if (k<words) for (unsigned i=0;i<3;++i) for (unsigned j=0;j<3;++j)
        if (!Complete || (i<2 && j<2)) totals[3*i+j]+=population_count<true>(a[i*words+k]&b[j*words+k]);
    Counts result;
    for (unsigned i=0;i<3;++i) for (unsigned j=0;j<3;++j) {
        const double value=static_cast<double>(totals[3*i+j]);
        result.cells[3*i+j]=value;
        if (!Complete) { result.rows[i]+=value; result.cols[j]+=value; result.n+=value; }
    }
    return result;
}
#endif
#if defined(GWLD_X86_SIMD)
// AVX2 has no vector POPCNT: count low/high nibbles through a byte lookup.
// Flush byte accumulators after 31 vectors (31*8=248), before they wrap.
template<bool Complete=false>
__attribute__((target("avx2")))
GWLD_FORCE_INLINE Counts count_packed_avx2(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    const __m256i lookup=_mm256_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
                                         0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
    const __m256i mask=_mm256_set1_epi8(15),zero=_mm256_setzero_si256();
    Counts result;
    for (unsigned i=0;i<(Complete ? 2u:3u);++i) for (unsigned j=0;j<(Complete ? 2u:3u);++j) {
        __m256i sum=_mm256_setzero_si256(); std::size_t k=0;
        while (words-k>=4) {
            const auto end=k+4*std::min<std::size_t>((words-k)/4,31);
            __m256i bytes=_mm256_setzero_si256();
            for (;k<end;k+=4) {
                const auto x=_mm256_and_si256(
                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a+i*words+k)),
                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b+j*words+k)));
                const auto low=_mm256_shuffle_epi8(lookup,_mm256_and_si256(x,mask));
                const auto high=_mm256_shuffle_epi8(lookup,_mm256_and_si256(_mm256_srli_epi16(x,4),mask));
                bytes=_mm256_add_epi8(bytes,_mm256_add_epi8(low,high));
            }
            sum=_mm256_add_epi64(sum,_mm256_sad_epu8(bytes,zero));
        }
        std::uint64_t lanes[4]; _mm256_storeu_si256(reinterpret_cast<__m256i*>(lanes),sum);
        std::uint64_t total=lanes[0]+lanes[1]+lanes[2]+lanes[3];
        for (;k<words;++k) total+=population_count<false>(a[i*words+k]&b[j*words+k]);
        const double value=static_cast<double>(total);
        result.cells[3*i+j]=value;
        if (!Complete) { result.rows[i]+=value; result.cols[j]+=value; result.n+=value; }
    }
    return result;
}
#endif
#if defined(GWLD_X86_AVX512)
// AVX512F alone does not provide vector population counts; require VPOPCNTDQ.
template<bool Complete=false>
__attribute__((target("avx512f,avx512vpopcntdq")))
GWLD_FORCE_INLINE Counts count_packed_avx512(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    Counts result;
    for (unsigned i=0;i<(Complete ? 2u:3u);++i) for (unsigned j=0;j<(Complete ? 2u:3u);++j) {
        __m512i sum=_mm512_setzero_si512(); std::size_t k=0;
        for (;words-k>=8;k+=8) {
            const auto x=_mm512_and_si512(_mm512_loadu_si512(a+i*words+k),_mm512_loadu_si512(b+j*words+k));
            sum=_mm512_add_epi64(sum,_mm512_popcnt_epi64(x));
        }
        std::uint64_t lanes[8]; _mm512_storeu_si512(lanes,sum);
        std::uint64_t total=0; for (auto lane:lanes) total+=lane;
        for (;k<words;++k) total+=population_count<false>(a[i*words+k]&b[j*words+k]);
        const double value=static_cast<double>(total);
        result.cells[3*i+j]=value;
        if (!Complete) { result.rows[i]+=value; result.cols[j]+=value; result.n+=value; }
    }
    return result;
}
#endif

// Four doubles use the homogeneous floating-point return ABI on ARM64.
// Keep large Counts objects out of the complete-pair kernel call boundary.
struct FourCounts { double c00,c01,c10,c11; };
GWLD_FORCE_INLINE FourCounts four_counts(const Counts& c) {
    return {c.cells[0],c.cells[1],c.cells[3],c.cells[4]};
}
typedef FourCounts (*FourCountKernel)(const std::uint64_t*,const std::uint64_t*,std::size_t);
typedef Counts (*CountKernel)(const std::uint64_t*,const std::uint64_t*,std::size_t);
inline CountKernel select_count_kernel(std::size_t words=0) {
    (void)words;
    #if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
    #if defined(GWLD_X86_AVX512)
    if (words>=8 && __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vpopcntdq"))
        return count_packed_avx512<>;
    #endif
    if (words>=16 && __builtin_cpu_supports("avx2")) return count_packed_avx2<>;
    return __builtin_cpu_supports("popcnt") ? count_packed_native : count_packed_portable;
    #elif defined(_MSC_VER) && defined(_M_X64)
    int registers[4]; __cpuid(registers,1);
    return (registers[2]&(1<<23)) ? count_packed_native : count_packed_portable;
    #elif defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    return count_packed_neon<>;
    #elif defined(__GNUC__) || defined(__clang__)
    return count_packed_native; // compiler lowers the builtin for the target ISA
    #else
    return count_packed_portable;
    #endif
}

inline FourCounts count_four_portable(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    return four_counts(packed_count_kernel<false,true>(a,b,words));
}
#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
__attribute__((target("popcnt")))
#endif
inline FourCounts count_four_native(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    return four_counts(packed_count_kernel<true,true>(a,b,words));
}
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
inline FourCounts count_four_neon(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    return four_counts(count_packed_neon<true>(a,b,words));
}
#endif
#if defined(GWLD_X86_SIMD)
__attribute__((target("avx2")))
inline FourCounts count_four_avx2(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    return four_counts(count_packed_avx2<true>(a,b,words));
}
#endif
#if defined(GWLD_X86_AVX512)
__attribute__((target("avx512f,avx512vpopcntdq")))
inline FourCounts count_four_avx512(const std::uint64_t* a,const std::uint64_t* b,std::size_t words) {
    return four_counts(count_packed_avx512<true>(a,b,words));
}
#endif
inline FourCountKernel select_complete_count_kernel(CountKernel selected) {
    if (selected==count_packed_portable) return count_four_portable;
    #if defined(GWLD_X86_SIMD)
    if (selected==count_packed_avx2<>) return count_four_avx2;
    #endif
    #if defined(GWLD_X86_AVX512)
    if (selected==count_packed_avx512<>) return count_four_avx512;
    #endif
    if (selected==count_packed_native) return count_four_native;
    #if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    return count_four_neon;
    #else
    return count_four_native;
    #endif
}

class PackedGenotypes {
    std::size_t words=0;
    std::vector<std::uint64_t> bits;
    CountKernel kernel=count_packed_portable;
    FourCountKernel complete_kernel=count_four_portable;
public:
    template<class Matrix>
    void build(const Matrix& matrix,CountKernel selected=nullptr) {
        words=matrix.n_rows/64+(matrix.n_rows%64!=0);
        if (!words || !matrix.n_cols) { bits.clear(); return; }
        if (words>bits.max_size()/3 || matrix.n_cols>bits.max_size()/(3*words))
            throw std::length_error("Packed genotype matrix is too large");
        bits.assign(matrix.n_cols*3*words,0); kernel=selected ? selected : select_count_kernel(words);
        complete_kernel=select_complete_count_kernel(kernel);
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
    Counts count_complete(std::size_t i,std::size_t j,const Summary& a,const Summary& b) const {
        const FourCounts cells=complete_kernel(bits.data()+i*3*words,bits.data()+j*3*words,words);
        Counts result;
        result.cells[0]=cells.c00; result.cells[1]=cells.c01;
        result.cells[3]=cells.c10; result.cells[4]=cells.c11;
        // With no missing samples, four intersections determine the entire 3x3 table.
        for (unsigned k=0;k<3;++k) { result.rows[k]=a.counts[k]; result.cols[k]=b.counts[k]; }
        result.n=a.total;
        result.cells[2]=a.counts[0]-result.cells[0]-result.cells[1];
        result.cells[5]=a.counts[1]-result.cells[3]-result.cells[4];
        result.cells[6]=b.counts[0]-result.cells[0]-result.cells[3];
        result.cells[7]=b.counts[1]-result.cells[1]-result.cells[4];
        result.cells[8]=a.counts[2]-result.cells[6]-result.cells[7];
        return result;
    }
    Counts count(std::size_t i,std::size_t j) const {
        return kernel(bits.data()+i*3*words,bits.data()+j*3*words,words);
    }
};
} // namespace GWLD
#undef GWLD_FORCE_INLINE
