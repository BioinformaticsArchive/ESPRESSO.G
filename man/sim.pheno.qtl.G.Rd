\name{sim.pheno.qtl.G}
\alias{sim.pheno.qtl.G}
\title{Simulates continuous outcome data}
\usage{
sim.pheno.qtl.G(num.subjects = 10000, pheno.mean = 0, pheno.sd = 1,
  genotype = NULL, geno.efkt = 0.25)
}
\arguments{
  \item{num.subjects}{number of subjects to generate.}

  \item{pheno.mean}{statistical mean}

  \item{pheno.sd}{standard deviation}

  \item{genotype}{a vector that represents the exposure
  data}

  \item{geno.efkt}{effect size of related to the 'at risk'
  allele.}
}
\value{
a binary vector that represents the phenotype data.
}
\description{
The function uses the effects data of the genetic
determinant to construct a linear predictor(LP). The
outcome is a normally distributed variable generated with a
mean equal to the LP and a standard deviation of 1.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

