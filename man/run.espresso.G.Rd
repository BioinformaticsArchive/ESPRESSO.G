\name{run.espresso.G}
\alias{run.espresso.G}
\title{Runs a full ESPRESSO analysis}
\usage{
run.espresso.G(simulation.params = NULL, pheno.params = NULL,
  geno.params = NULL, scenarios2run = 1)
}
\arguments{
  \item{simulation.params}{general parameters for the
  scenario(s) to analyse}

  \item{pheno.params}{paramaters for the outcome variables}

  \item{geno.params}{parameters for the genetic
  determinant}

  \item{scenarios2run}{the indices of the scenarios one
  wish to analyse}
}
\value{
a summary table that contains both the input parameters and
the results of the analysis
}
\description{
This function calls the functions required to run a full
ESPRESSO analysis where the model consists of an outcome
(binary or continuous) determinant by a binary or
continuous environmental determinant.
}
\details{

}
\examples{
{

# load the table that hold the input parameters; each of the table hold parameters for 4 scenarios:
# scenario 1: a binary outcome determined by a binary exposure
# scenario 2: a binary outcome determined by a quantitative and normally distributed exposure
# scenario 3: a quantitative outcome determined by a binary exposure
# scenario 4: a quantitative outcome determined by a quantitative and normally distributed exposure
data(simulation.params)
data(pheno.params)
data(geno.params)

# run the function for the first two scenarios (by default only the first scenario is ran)
run.espresso.G(simulation.params, pheno.params, geno.params, scenarios2run=c(1,2))

# run the function for the last two scenarios (by default only the first scenario is ran)
run.espresso.G(simulation.params, pheno.params, geno.params, scenarios2run=c(3,4))
}
}
\author{
Amadou Gaye
}

