\name{extract_mod_metadata}
\alias{extract_mod_metadata}
\title{Extract \code{@dynhr:} metadata blocks from a .mod file}
\usage{
extract_mod_metadata(mod_file_or_lines)
}
\arguments{
\item{mod_file_or_lines}{Path to a Dynare \code{.mod} file, or a character
  vector of lines already read from one.}
}
\value{
Named list of metadata blocks parsed from \code{@dynhr:*} comment
  sections in the \code{.mod} file:  \code{shock_categories},
  \code{benchmarks}, \code{narratives}, \code{variable_labels}, etc.
  Returns an empty list if no metadata blocks are found.
}
\description{
Reads a \code{.mod} file and extracts the structured \code{@dynhr:}
metadata blocks that extend the Dynare format.  These blocks drive the
diagnostic battery (D-suite) without requiring a separate configuration
file.
}
\seealso{
\code{\link{parse_mod}}
}
