#' 
#' @title Adds some error to genotype data
#' @description Simulates errors and adds it to the true data to obtain observed data. 
#'  The alleles simulated by the function sim.geno.data are randomly misclassified and used to 
#'  form new genotypes that represent the observed genotypes.
#' @param allele.A Allele A
#' @param allele.B Allele B
#' @param geno.model genetic model; binary=0 and additive=1
#' @param MAF minor allele frequency of the SNP (in ESPRESSO this is the frequency of the 'at risk' allele)
#' @param geno.error a vector with two values, the misclassififcation rates related to the sensitivity and 
#' specificity of the assessment of the alleles, i.e. 1-sensitivity and 1-specificity
#' @return a dataframe that contains the below data:
#' \code{observed.genotype} observed genotypes
#' \code{observed.allele.A} observed A alleles
#' \code{observed.allele.B} observed B alleles
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # generate data for a binary SNP with a MAF of 0.1
#' geno.elts <- sim.geno.data(num.obs=10000, geno.model=0, MAF=0.1)
#' allele.A <- geno.elts$allele.A
#' allele.B <- geno.elts$allele.B
#' 
#' # generate the observed genotypes by inducing some misclassification determined 
#' # by the sensitivity and specificity of the assessment of the alleles
#' observed.geno.elts <- get.obs.geno(allele.A, allele.B, geno.model=0, MAF=0.1, geno.error=c(0.05, 0.05))
#' observed.geno.data <- observed.geno.elts$observed.genotype
#' 
#' }
#' 
get.obs.geno <- function (allele.A=NULL, allele.B=NULL, geno.model=0, MAF=0.1, geno.error=c(0.05, 0.05)){
  
  # IF ALLELE DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
  if(is.null(allele.A)){
    cat("\n\n ALERT!\n")
    cat(" No allele data found for the first allele.\n")
    cat(" Check the argument 'allele.A'\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  if(is.null(allele.B)){
    cat("\n\n ALERT!\n")
    cat(" No allele data found for the second allele.\n")
    cat(" Check the argument 'allele.B'\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
			
  mean.add <- (2*MAF*(1-MAF)+2*(MAF^2))
  mean.bin <- (2*MAF*(1-MAF)+(MAF^2))
  true.allele.A <- allele.A
  true.allele.B <- allele.B
  
  observed.allele.A <- misclassify(binary.vector=true.allele.A, error.1.0=geno.error[1], error.0.1=geno.error[2])
  observed.allele.B <- misclassify(binary.vector=true.allele.B, error.1.0=geno.error[1], error.0.1=geno.error[2])
  observed.genotype <- observed.allele.A+observed.allele.B
  
  if(geno.model==0){
    genotyp.U <- observed.genotype > 0
    observed.genotype <- genotyp.U - mean.bin
  }
  if(geno.model==1){
    genotyp.U <- observed.genotype
    observed.genotype <- genotyp.U - mean.add
  }
  
  df <- data.frame(observed.genotype, observed.allele.A, observed.allele.B)
}

