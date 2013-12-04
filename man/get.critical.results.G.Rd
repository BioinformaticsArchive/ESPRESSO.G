\name{get.critical.results.G}
\alias{get.critical.results.G}
\title{Provides a summary of the main results}
\usage{
  get.critical.results.G(scenario = 1, pheno.model = 0,
    geno.model = 1, sample.sizes.required = NULL,
    empirical.power = 0.8, modelled.power = 0.8,
    mean.beta = NULL)
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
  a list that contains the following items:
  \code{environment.model} model of the environmental
  determinant \code{number.of.cases.required} number of
  cases required to achieve the desired power under a
  binary outcome model. \code{number.of.controls.required}
  number of controls required to achieve the desired power
  under a binary outcome model.
  \code{number.of.subjects.required} number of subjects
  required to achieve the desired power under a
  quantatative outcome model. \code{empirical.power}
  estimated empirical power under each model.
  \code{modelled.power} power achieved under each model
  with specified sample size. \code{estimated.OR} esitmated
  odds-ratios due to shrinkage toward the null resulting
  from misclassification.
}
\description{
  Gets the number of cases and controls or subjects and the
  empirical and theoretical power under each model and
  prints a summary on the screen.
}
\examples{
{

# This example assumes results obtained after running an analysis with a quantitative outcome
# determined by a binary SNP

# scenario number
j <- 1

# outcome model
pheno.model <- 1

# model of the environmental exposure
geno.model <- 0

# Estimated number of subjects required to achieve 80\% power
sample.sizes.required <- 10000

# estimated values for empirical power from the input sample size
empirical.power <- 0.2

# estimated values for modelled power from the input sample size
modelled.power <- 0.19

# mean beta value for the single environmenatl determinant
mean.beta <- 0.20

# return critical results and print a summary
res <- get.critical.results.G(j, pheno.model, geno.model, sample.sizes.required,
empirical.power, modelled.power, mean.beta)
}
}
\author{
  Amadou Gaye
}
