\name{sim.pheno.qtl.G}
\alias{sim.pheno.qtl.G}
\title{Simulates continuous outcome data}
\usage{
sim.pheno.qtl.G(num.subjects = 10000, genotype = NULL, geno.efkt = 0.25)
}
\arguments{
  \item{num.subjects}{number of subjects to generate.}

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
determinant to construct a linear predictor(LP).  The
outcome is a normally distributed variable generated with a
mean equal to the LP and a standard deviation of 1.
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

# Generate quantitative phenotype statuses using the exposure data obtained above and an effect size of 0.25
# for the 'at risk' genotype
pheno.data <- sim.pheno.qtl.G(num.subjects=10000, genotype=geno.data, geno.efkt=0.25)
}
}
\author{
Amadou Gaye
}

