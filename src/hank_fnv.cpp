// A 64-bit FNV-1a hash of a raw vector.
//
// Exists so hank_het3_fingerprint() can hash a serialised block WITHOUT a new
// dependency: `digest` is only in Suggests, and a public identity function
// that silently changes behaviour depending on whether an optional package is
// installed is worse than no function at all.
//
// FNV-1a is not cryptographic and is not claimed to be. The job here is to
// tell APART two households that differ somewhere in a few million doubles --
// accidental, not adversarial, collisions. For that it is well past
// sufficient, and it is deterministic across platforms because the input is
// R's own XDR (big-endian) serialisation, so a fingerprint computed on one
// machine matches one computed on another.

#include <Rcpp.h>
#include <cstdint>
#include <cstdio>
using namespace Rcpp;

// [[Rcpp::export]]
std::string hank_fnv1a64_cpp(RawVector x) {
  uint64_t h = 14695981039346656037ULL;          // FNV offset basis
  const uint64_t prime = 1099511628211ULL;       // FNV prime
  const R_xlen_t n = x.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    h ^= static_cast<uint64_t>(x[i]);
    h *= prime;
  }
  char buf[17];
  std::snprintf(buf, sizeof(buf), "%016llx",
                static_cast<unsigned long long>(h));
  return std::string(buf);
}
