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
\examples{
{

# Generate data for a binary SNP with a MAF of 0.1
geno.elts <- sim.geno.data(num.obs=10000, geno.model=0, MAF=0.1)
geno.data <- geno.elts$genotype

# simulate subject effect data
s.effect <- sim.subject.data(num.obs=10000, baseline.OR=10)

# Generate binary phenotype statuses using the exposure data obtained above and an OR of 1.5 for the 'at risk' genotype
# and a prevalence of 0.1 for the binary outcome
pheno.data <- sim.pheno.bin.G(num.obs=10000,disease.prev=0.1, genotype=geno.data, subject.effect.data=s.effect, geno.OR=1.5)

}
}
\author{
Amadou Gaye
}

