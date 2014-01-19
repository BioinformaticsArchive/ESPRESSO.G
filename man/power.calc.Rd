\name{power.calc}
\alias{power.calc}
\title{Calculates the empirical and theoretical power}
\usage{
power.calc(pval = 1e-04, z.values = NULL, mean.model.z = NULL)
}
\arguments{
  \item{pval}{cut-off p-value defining statistical
  significance.}

  \item{z.values}{z-statistic of the determinant.}

  \item{mean.model.z}{mean z-statistic of the environmental
  determinant.}
}
\value{
a list that contains the computed empirical power and
theoretical power.
}
\description{
The function determines the empirical and theoretical
power. The empirical power is the proportion of simulations
in which the z-statistic for the parameter of interest
exceeds the z-statistic for the desured level if
statistical significance. The theoretical power is the
power of the study.
}
\details{

}
\examples{
{

# generate the arguments to pass on to the function (typically the output of the function glm.analysis.G)
# load a table of quantitative outcome and binary SNP data
data(observed.data.G)

# run a logistic regression (outcome is binary)
glm.estimates <- glm.analysis.G(pheno.model=1, observed.data.G)
beta <- glm.estimates$beta
se <- glm.estimates$se
modz <- beta/se
zvals <-  glm.estimates$z

# use the function; in a typical ESPRESSO process 'modz' would be the ratio of the average beta over the average
# se across the stored beta and se values (after each run the estimates of the glm analysis are stored) and 'zvals'
# would be the vector of z values obtained after the iterations.
power.values <- power.calc(z.values=zvals, mean.model.z=modz)
}
}
\author{
Amadou Gaye
}

