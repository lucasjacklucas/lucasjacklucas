\name{latentvar}
\alias{latentvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple summary of OLS model results
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
This function runs an OLS model and spits out useful information: a results table, a model plot, and a summary of the variables in the model. It also saves base and tidy versions of the model, as well as the model plot, in a list.
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

