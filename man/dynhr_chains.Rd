\name{dynhr_chains}
\alias{dynhr_chains}
\docType{class}
\title{S3 class for posterior sample chains}
\description{
An S3 class returned by \code{\link{dynhr_mcmc}}, \code{\link{smc}}, and
\code{\link{nuts}}.  Wraps the raw sampler output with consistent slots
and \code{print} / \code{summary} methods.
}
\section{Slots}{
\describe{
  \item{\code{chain}}{Numeric matrix (\eqn{n_{\text{draws}} \times n_p})
    of post-warmup posterior draws.  Column names are parameter names.}
  \item{\code{sampler}}{Character: one of \code{"rwmh"}, \code{"smc"}, or
    \code{"nuts"}.}
  \item{\code{acceptance_rate}}{Mean acceptance rate: the Metropolis
    acceptance rate for RWMH; for NUTS the mean post-warmup NUTS acceptance
    statistic (the quantity \code{target_accept} targets; each chain also
    carries \code{accept_stat} and \code{move_rate}, the fraction of
    iterations whose draw moved).}
  \item{\code{log_marginal_lik}}{Log marginal likelihood estimate (SMC
    only).}
  \item{\code{n_divergent}}{Number of divergent transitions (NUTS only).}
  \item{\code{elapsed_secs}}{Wall-clock seconds.}
}
}
\seealso{
\code{\link{dynhr_mcmc}}, \code{\link{smc}}, \code{\link{nuts}}
}
