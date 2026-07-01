\name{dynhr_chains}
\alias{dynhr_chains}
\docType{class}
\title{S3 class for posterior sample chains}
\description{
An S3 class returned by \code{\link{mcmc}}, \code{\link{smc}}, and
\code{\link{nuts}}.  Wraps the raw sampler output with consistent slots
and \code{print} / \code{summary} methods.
}
\section{Slots}{
\describe{
  \item{\code{chain}}{Numeric matrix (\eqn{n_{\text{draws}} \times n_p})
    of post-warmup posterior draws.  Column names are parameter names.}
  \item{\code{sampler}}{Character: one of \code{"rwmh"}, \code{"smc"}, or
    \code{"nuts"}.}
  \item{\code{acceptance_rate}}{Mean acceptance rate (RWMH only).}
  \item{\code{log_marginal_lik}}{Log marginal likelihood estimate (SMC
    only).}
  \item{\code{n_divergent}}{Number of divergent transitions (NUTS only).}
  \item{\code{elapsed_secs}}{Wall-clock seconds.}
}
}
\seealso{
\code{\link{mcmc}}, \code{\link{smc}}, \code{\link{nuts}}
}
