\name{sim.QTL.data.G}
\alias{sim.QTL.data.G}
\title{Simulates subjects for a continuous outcome}
\usage{
sim.QTL.data.G(numsubjects = 10000, ph.mean = 0, ph.sd = 1,
  geno.model = 0, MAF = 0.1, geno.efkt = 0.25, pheno.reliability = 0.9)
}
\arguments{
  \item{numsubjects}{number of subjects to simulate.}

  \item{ph.mean}{statistical mean}

  \item{ph.sd}{standard deviation}

  \item{geno.model}{genetic model; binary=0 and
  additive=1.}

  \item{MAF}{minor allele frequencies of the genetic
  variant.}

  \item{geno.efkt}{effect size of the 'at risk' genotype.}

  \item{pheno.reliability}{reliability of the assessment
  for a quantitative outcome.}
}
\value{
a matrix that holds the outcome (\code{phenotype}) and
exposure (\code{genotype}) data and the SNP alleles used to
construct the genotypes.
}
\description{
Generates the specified number of subjects for a
quantitative outcome.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

