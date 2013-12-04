\name{glm.analysis.G}
\alias{glm.analysis.G}
\title{Fits a Generalized Linear Model}
\usage{
  glm.analysis.G(pheno.model = 0, observed.data = NULL)
}
\arguments{
  \item{pheno.model}{outcome type: binary=0 and normal=1.}

  \item{observed.data}{a dataframe that contains the
  observerd outcome and covariate data.}
}
\value{
  a vector that holds the beta, standard-error and
  z-statistic of each of the covariates.
}
\description{
  This function fits a conventional unconditional logistic
  regression model on the observed outcome and exposure
  data. In a Typical ESPRESSO the function is called at
  each iteration and the estimates stored.
}
\examples{
{

# load a table that contains observed quantitative outcome and binary SNP data
data(observed.data.G)

# run a logistic regression (outcome is binary)
glm.estimates <- glm.analysis.G(pheno.model=1, observed.data.G)

}
}
\author{
  Amadou Gaye
}

