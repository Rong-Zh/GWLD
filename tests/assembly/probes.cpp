// Stable names expose optimized kernels for inspection without Armadillo/R.
#include "gwld/packed.hpp"
extern "C" {
__attribute__((noinline)) void gwld_constant_probe(const GWLD::Summary* a,const GWLD::Summary* b,GWLD::Counts* out) {
    *out=GWLD::count_constant(*a,*b);
}
__attribute__((noinline)) void gwld_pack_probe(const int* g,std::size_t rows,std::size_t words,std::uint64_t* out) {
    GWLD::packing_detail::column(g,rows,words,out);
}
__attribute__((noinline)) void gwld_summary_probe(const int* g,std::size_t n,GWLD::Summary* out) {
    *out=GWLD::summarize(g,n);
}
__attribute__((noinline)) void gwld_scalar_probe(const int* a,const int* b,std::size_t n,GWLD::Counts* out) {
    *out=GWLD::count(a,b,n);
}
__attribute__((noinline)) GWLD::CountKernel gwld_dispatch_probe(std::size_t words) {
    return GWLD::select_count_kernel(words);
}
__attribute__((noinline)) GWLD::FourCountKernel gwld_complete_dispatch_probe(GWLD::CountKernel p) {
    return GWLD::select_complete_count_kernel(p);
}
#if defined(__aarch64__)
__attribute__((noinline)) GWLD::FourCounts gwld_neon_four_probe(const std::uint64_t* a,const std::uint64_t* b,std::size_t n) {
    return GWLD::count_four_neon(a,b,n);
}
#endif
__attribute__((noinline)) double gwld_likelihood_probe(double x,double p,double q,const GWLD::Counts* counts) {
    return GWLD::nll(x,p,q,*counts);
}
}
