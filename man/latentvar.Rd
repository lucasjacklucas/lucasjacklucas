\name{latentvar}
\alias{latentvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian factor analysis for latent variable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
latentvar(data.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.name}{
%%     ~~Describe \code{data.name} here~~
}
}
\details{
This function estimates a latent variable using a set of indicator variables. Note that the first column in the input data frame should be an identifier code (such as a respondent ID), and all other columns should be indicator variables for the latent measure.
}
\value{
\item{results }{Data frame containing latent estimate for each observation, including latent estimate as well as lower and upper credible intervals.}
\item{betas }{Discrimination parameters for each item.}
\item{beta.plot }{Simple plot of discrimination parameters.}
\item{alldraws }{Complete posterior draws for each observation/respondent.}
\item{convergence }{Convergence information for model}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

