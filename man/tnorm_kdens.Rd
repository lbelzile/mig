% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bandwidth.R
\name{tnorm_kdens}
\alias{tnorm_kdens}
\title{Truncated Gaussian kernel density estimator}
\usage{
tnorm_kdens(x, newdata, Sigma, beta, log = TRUE, ...)
}
\arguments{
\item{x}{\code{n} by \code{d} matrix of quantiles}

\item{newdata}{matrix of new observations at which to evaluated the kernel density}

\item{Sigma}{scale matrix}

\item{beta}{\code{d} vector \eqn{\boldsymbol{\beta}} defining the half-space through \eqn{\boldsymbol{\beta}^{\top}\boldsymbol{\xi}>0}}

\item{log}{logical; if \code{TRUE}, returns log probabilities}

\item{...}{additional arguments, currently ignored}
}
\value{
a vector containing the value of the kernel density at each of the \code{newdata} points
}
\description{
Given a data matrix over a half-space defined by \code{beta},
compute the log density of the asymmetric truncated Gaussian kernel density estimator,
taking in turn an observation as location vector.
}
