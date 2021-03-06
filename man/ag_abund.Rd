% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agg_funcs.R
\name{ag_abund}
\alias{ag_abund}
\title{Aggregate abundances based on a factor variable}
\usage{
ag_abund(x, agg.var)
}
\arguments{
\item{x}{A data set containing columns '\code{agg.var}', '\code{y_pred}'
and '\code{y_real}'. The last two are obtained from a calls to
\code{agTrendNimble::fit_ssl_nimble}.}

\item{agg.var}{Variable used for aggregation}
}
\description{
Aggregate abundances based on a factor variable
}
