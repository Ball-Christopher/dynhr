// What the compiler actually did for dynhr's own shared library.
//
// R's `CMD config CXXFLAGS` reports R's DEFAULT flags, which say nothing
// about how this DLL was built when the user installed with a personal
// Makevars (or a distribution that overrides flags). The predefined macros
// below are evaluated at compile time of THIS translation unit, so they
// describe the flags the package's C++ sources were really compiled with.

#include <Rcpp.h>
#include <string>

//' Compile-time facts about dynhr's own shared library
//'
//' Returns the predefined compiler macros (optimisation level, fast-math,
//' FMA / SIMD instruction-set flags, compiler identity, C++ standard) as seen
//' when the package's C++ was compiled. Internal: surfaced to users through
//' \code{\link{dynhr_system_info}}.
//'
//' @return A named list.
//' @keywords internal
//' @noRd
// [[Rcpp::export(name = ".dynhr_build_info")]]
Rcpp::List dynhr_build_info_cpp() {
  using Rcpp::Named;
#ifdef __OPTIMIZE__
  const bool optimize = true;
#else
  const bool optimize = false;
#endif
#ifdef __OPTIMIZE_SIZE__
  const bool optimize_size = true;
#else
  const bool optimize_size = false;
#endif
#ifdef __FAST_MATH__
  const bool fast_math = true;
#else
  const bool fast_math = false;
#endif
#ifdef __FMA__
  const bool fma = true;
#else
  const bool fma = false;
#endif
#ifdef __AVX__
  const bool avx = true;
#else
  const bool avx = false;
#endif
#ifdef __AVX2__
  const bool avx2 = true;
#else
  const bool avx2 = false;
#endif
#ifdef __AVX512F__
  const bool avx512f = true;
#else
  const bool avx512f = false;
#endif
#ifdef __SSE2__
  const bool sse2 = true;
#else
  const bool sse2 = false;
#endif
#ifdef __ARM_NEON
  const bool neon = true;
#else
  const bool neon = false;
#endif
#ifdef __aarch64__
  const bool aarch64 = true;
#else
  const bool aarch64 = false;
#endif
#ifdef FP_FAST_FMA
  const bool fp_fast_fma = true;
#else
  const bool fp_fast_fma = false;
#endif

  std::string compiler = "unknown";
  std::string compiler_version = "";
#if defined(__clang__)
  compiler = "clang";
#  ifdef __clang_version__
  compiler_version = __clang_version__;
#  endif
#elif defined(__GNUC__)
  compiler = "gcc";
  compiler_version = std::to_string(__GNUC__) + "." +
    std::to_string(__GNUC_MINOR__) + "." + std::to_string(__GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER)
  compiler = "msvc";
  compiler_version = std::to_string(_MSC_VER);
#endif

  return Rcpp::List::create(
    Named("optimize") = optimize,
    Named("optimize_size") = optimize_size,
    Named("fast_math") = fast_math,
    Named("fma") = fma,
    Named("fp_fast_fma") = fp_fast_fma,
    Named("avx") = avx,
    Named("avx2") = avx2,
    Named("avx512f") = avx512f,
    Named("sse2") = sse2,
    Named("neon") = neon,
    Named("aarch64") = aarch64,
    Named("compiler") = compiler,
    Named("compiler_version") = compiler_version,
    Named("cplusplus") = static_cast<double>(__cplusplus));
}
