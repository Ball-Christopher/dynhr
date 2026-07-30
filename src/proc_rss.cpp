// Peak resident-set-size reporter, for the three-asset run manifest (A8).
//
// This is a PEAK-SINCE-PROCESS-START figure -- monotone non-decreasing for
// the life of the R session -- not "the memory this call used". It answers
// "how much RAM has this process touched so far", which is what a release
// manifest's memory-budget check needs.
//
// Units are the classic cross-platform trap here: POSIX getrusage() reports
// ru_maxrss in BYTES on macOS/Darwin but KILOBYTES on Linux. Get the branch
// wrong and every number is off by exactly 1024x -- silently, because both
// look like plausible RSS figures in isolation. Windows has no getrusage();
// GetProcessMemoryInfo()'s PeakWorkingSetSize is already in bytes.

#include <Rcpp.h>

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__unix__) || defined(__unix) || defined(__APPLE__)
#include <sys/resource.h>
#endif

//' Peak resident set size of the current process, in bytes
//'
//' Reads the operating system's own peak-RSS counter for the calling R
//' process. This is a PEAK-SINCE-PROCESS-START figure -- monotone
//' non-decreasing across the whole R session -- not the memory consumed by
//' this call or by any single solve. It exists to feed the three-asset run
//' manifest's peak-memory column (see \code{\link{hank_het3_manifest}}).
//'
//' \strong{Platform coverage.} POSIX (macOS/Linux/BSD) uses
//' \code{getrusage(RUSAGE_SELF, \&ru)} and reads \code{ru_maxrss}. The units
//' of that field are NOT portable: macOS/Darwin reports bytes, Linux
//' reports kilobytes, and this function corrects for that internally so its
//' OWN return value is always bytes. Windows uses
//' \code{GetProcessMemoryInfo()}'s \code{PeakWorkingSetSize}, which is
//' already in bytes. On any platform where none of these APIs is available,
//' the function returns \code{NA_real_} rather than guess.
//'
//' @return A single \code{double}: peak RSS in bytes, or \code{NA_real_}
//'   where the platform cannot report it.
//' @seealso \code{\link{hank_het3_manifest}}
//' @keywords internal
// [[Rcpp::export(name=".hank_peak_rss")]]
double hank_peak_rss() {
#if defined(_WIN32)
  PROCESS_MEMORY_COUNTERS pmc;
  if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
    return static_cast<double>(pmc.PeakWorkingSetSize);
  }
  return NA_REAL;
#elif defined(__unix__) || defined(__unix) || defined(__APPLE__)
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0) {
    return NA_REAL;
  }
  double maxrss = static_cast<double>(ru.ru_maxrss);
#ifdef __APPLE__
  // Darwin: ru_maxrss is already in bytes.
  return maxrss;
#else
  // Linux (and most other unices): ru_maxrss is in kilobytes.
  return maxrss * 1024.0;
#endif
#else
  return NA_REAL;
#endif
}
