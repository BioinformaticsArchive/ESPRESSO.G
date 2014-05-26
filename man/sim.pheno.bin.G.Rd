\name{sim.pheno.bin.G}
\alias{sim.pheno.bin.G}
\title{Generates phenotype statuses}
\usage{
sim.pheno.bin.G(num.obs = 10000, disease.prev = 0.1, genotype = NULL,
  subject.effect.data = NULL, geno.OR = 1.5)
}
\arguments{
  \item{num.obs}{number of observations to generate per
  iteration.}

  \item{disease.prev}{prevalence of the binary outcome.}

  \item{genotype}{a vector that represents the exposure
  data.}

  \item{subject.effect.data}{subject effect data, reflects
  the heterogenity in baseline disease risk.}

  \item{geno.OR}{odds ratio related to the 'at risk'
  genotype.}
}
\value{
a binary vector that represents the phenotype data.
}
\description{
Generates affected and non-affected subjects using the
genotypes.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

