\name{sim.QTL.data.G}
\alias{sim.QTL.data.G}
\title{Simulates subjects for a continuous outcome}
\usage{
sim.QTL.data.G(numsubjects = 10000, geno.model = 0, MAF = 0.1,
  geno.efkt = 0.25, pheno.reliability = 0.9)
}
\arguments{
  \item{numsubjects}{number of subjects to simulate.}

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
\examples{
{

# Example 1: generate 5000 subjects for a binary genetic variant with a MAF of 0.1 and an effect size of 0.25;
# the outcome is measured with a reliability of 0.9.
sim.matrix <- sim.QTL.data.G(numsubjects=5000, geno.model=0, MAF=0.1, geno.efkt=0.25, pheno.reliability=0.9)

# Example 2: generate 5000 subjects for an additive genetic variant with a  MAF of 0.1 and an effect size of 0.25;
# the outcome is measured with a reliability of 0.9.
sim.matrix <- sim.QTL.data.G(numsubjects=5000, geno.model=1, MAF=0.1, geno.efkt=0.25, pheno.reliability=0.9)

}
}
\author{
Amadou Gaye
}

