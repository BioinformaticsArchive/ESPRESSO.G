\name{get.obs.geno}
\alias{get.obs.geno}
\title{Adds some error to genotype data}
\usage{
  get.obs.geno(allele.A = NULL, allele.B = NULL,
    geno.model = 0, MAF = 0.1, geno.error = c(0.05, 0.05))
}
\arguments{
  \item{allele.A}{Allele A}

  \item{allele.B}{Allele B}

  \item{geno.model}{genetic model; binary=0 and additive=1}

  \item{MAF}{minor allele frequency of the SNP (in ESPRESSO
  this is the frequency of the 'at risk' allele)}

  \item{geno.error}{a vector with two values, the
  misclassififcation rates related to the sensitivity and
  specificity of the assessment of the alleles, i.e.
  1-sensitivity and 1-specificity}
}
\value{
  a dataframe that contains the below data:
  \code{observed.genotype} observed genotypes
  \code{observed.allele.A} observed A alleles
  \code{observed.allele.B} observed B alleles
}
\description{
  Simulates errors and adds it to the true data to obtain
  observed data.  The alleles simulated by the function
  sim.geno.data are randomly misclassified and used to form
  new genotypes that represent the observed genotypes.
}
\examples{
{

# generate data for a binary SNP with a MAF of 0.1
geno.elts <- sim.geno.data(num.obs=10000, geno.model=0, MAF=0.1)
allele.A <- geno.elts$allele.A
allele.B <- geno.elts$allele.B

# generate the observed genotypes by inducing some misclassification determined
# by the sensitivity and specificity of the assessment of the alleles
observed.geno.elts <- get.obs.geno(allele.A, allele.B, geno.model=0, MAF=0.1, geno.error=c(0.05, 0.05))
observed.geno.data <- observed.geno.elts$observed.genotype

}
}
\author{
  Amadou Gaye
}

