\name{get.critical.results.G}
\alias{get.critical.results.G}
\title{Provides a summary of the main results}
\usage{
get.critical.results.G(scenario = 1, pheno.model = 0, geno.model = 1,
  sample.sizes.required = NULL, empirical.power = 0.8,
  modelled.power = 0.8, mean.beta = NULL)
}
\arguments{
  \item{scenario}{scenario number.}

  \item{pheno.model}{outcome type; binary=0 and
  quantitative=1.}

  \item{geno.model}{genetic model; binary=0 and
  additive=1.}

  \item{sample.sizes.required}{number of cases and controls
  or number of subjects required to achieve the desired
  power.}

  \item{empirical.power}{estimated empirical power.}

  \item{modelled.power}{calculated theoretical power.}

  \item{mean.beta}{mean beta value of each of the
  determinants.}
}
\value{
a list that contains the sample size required to achieve
the desired power and the
}
\description{
Gets the number of cases and controls or subjects and the
empirical and theoretical power under each model and prints
a summary on the screen.
}
\details{

}
\author{
Gaye A.
}
\keyword{and}
\keyword{estimated}
\keyword{internal}
\keyword{modelled}
\keyword{power.}
\keyword{theoretical}

