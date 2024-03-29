% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_stats.R
\name{control_gramian}
\alias{control_gramian}
\title{Controllability Gramian}
\usage{
control_gramian(A, B, t = NA)
}
\arguments{
\item{A}{\eqn{n x n} matrix.}

\item{B}{\eqn{n x m} matrix.}

\item{t}{Either NA for infinite time Gramian, or a positive non-zero integer. Defaults to NA.}
}
\value{
The infinite time or finite time controllability Gramian
}
\description{
Compute the (infinite time) controllability Gramian for the discrete linear time invariant system described by \eqn{x(t+1) = Ax(t) + Bu(t)}. 
The infinite time controllability Gramian is the solution to the discrete Lyapunov equation \eqn{AWA^\prime-W = -BB^\prime}, while the finite time Gramian for time \eqn{T} is 
\deqn{W_t = \sum_{t = 0}^T A^tBB^\prime(A^\prime)^t}
}
\examples{

A = matrix(c(0,-3,-2,2,-2,1,-1,2,-1), 3,3)
B = diag(3)

#Infinite time Gramian
W_inf = control_gramian(A, B)

#4 time Gramian
W_4 = control_gramian(A,B,4)
}
