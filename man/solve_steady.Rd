\name{solve_steady}
\alias{solve_steady}
\title{Solve for the deterministic steady state}
\usage{
solve_steady(compiled, params, y0 = NULL, endo_names = NULL,
  exo_names = NULL, exo_init = NULL, max_iter = 1000L, tol = 1e-10,
  verbose = FALSE)
}
\arguments{
\item{compiled}{dynhr_compiled from \code{\link{compile_model}}.}
\item{params}{Named list or vector of parameter values.}
\item{y0}{Optional starting guess (named numeric vector).  If \code{NULL}
  the function builds its own starting point.}
\item{endo_names}{Character vector of endogenous variable names.}
\item{exo_names}{Character vector of exogenous variable names.}
\item{exo_init}{Optional named numeric vector of exogenous variable values
  to hold fixed when evaluating the static residuals (i.e. the steady-state
  level of the exogenous variables, from an \code{initval} block).  Names
  must match the model's \code{varexo} declarations; unnamed entries are
  matched in declaration order.  When \code{NULL} (default) all exogenous
  variables are treated as zero, which is the Dynare convention for models
  whose steady state is defined at \eqn{\varepsilon = 0}.}
\item{max_iter}{Maximum Newton iterations (default 1000).}
\item{tol}{Convergence tolerance on the residual norm (default 1e-10).}
\item{verbose}{Print progress messages (default \code{FALSE}).}
}
\value{
An object of class \code{dynhr_steady} with slots:
\describe{
  \item{\code{values}}{Named numeric vector of steady-state values.}
  \item{\code{residuals}}{Static residuals at the solution.}
  \item{\code{converged}}{Logical: did Newton converge?}
  \item{\code{iterations}}{Number of iterations taken.}
}
}
\description{
Solves the static residual system \eqn{F(\bar{y}) = 0} for the
deterministic steady state using Newton--Raphson with a backtracking
line search.  If the model has an \code{analytical_steady_state}
block that block is used instead.
}
\seealso{
\code{\link{compile_model}}, \code{\link{solve_perturbation}}
}
