// mvn_cdf.h -- the C++ accurate-mode log-CDF dispatch shared by mvn_cdf.cpp
// and pskf_filter.cpp.
#ifndef DYNHR_MVN_CDF_H
#define DYNHR_MVN_CDF_H

#include <vector>

// log Phi_q(x; 0, S) exactly as logcdf_ME_r(check = TRUE) evaluates it when
// every stage stays inside the C++ evaluators. x has q entries, S is q x q
// column-major. Returns false (out untouched) when the R code would leave
// them: a non-finite input, a bivariate term below p = 1e-3, a numerically
// singular lattice call, or a block larger than miwa_qmax.
bool logcdf_dispatch(const std::vector<double>& x, const std::vector<double>& S,
                     int q, double miwa_qmax, double& out);

#endif
