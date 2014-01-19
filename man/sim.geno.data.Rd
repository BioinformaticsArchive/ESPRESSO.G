\name{sim.geno.data}
\alias{sim.geno.data}
\title{Generates genotypes for a genetic variant}
\usage{
sim.geno.data(num.obs = 10000, geno.model = 0, MAF = 0.1)
}
\arguments{
  \item{num.obs}{number of observations to generate.}

  \item{geno.model}{genetic model; binary=0 and
  additive=1.}

  \item{MAF}{minor allele frequency of the SNP (in ESPRESSO
  this is the frequency of the 'at risk' allele)}
}
\value{
a dataframe that contains the following variables:
\code{allele.A} major allele \code{allele.B} minor allele
\code{genotype} genotype
}
\description{
Generates two alleles and combines them to form the
genotype of a SNP under a binary or additive genetic model.
}
\details{

}
\examples{
{

# Example 1: binary SNP
# Generate data for a binary SNP with a MAF of 0.1
geno.elts <- sim.geno.data(num.obs=10000, geno.model=0, MAF=0.1)
geno.data <- geno.elts$genotype

# Example 2: additive SNP
# Generate data for an additive SNP with a MAF of 0.1
geno.elts <- sim.geno.data(num.obs=10000, geno.model=1, MAF=0.1)
geno.data <- geno.elts$genotype

}
}
\author{
Amadou Gaye
}

