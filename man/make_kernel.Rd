% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{make_kernel}
\alias{make_kernel}
\title{Generate a Polynomial Kernel}
\usage{
make_kernel(degree, R)
}
\arguments{
\item{degree}{the degree of the polynomial of the kernel. Must be 4 or greater or NULL for
the uniform kernel. Note: (degree - 2) is the order of the kernel, meaning first non-zero
moment is x^(degree - 2).}

\item{R}{a positive number giving the support, (-R,R) of the kernel which is symmetric
and centered at 0.}
}
\description{
Generate a Polynomial Kernel
}
