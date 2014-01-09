#' 
#' @title Generates phenotype statuses
#' @description Generates affected and non-affected subjects using the genotypes.
#' @param num.obs number of observations to generate per iteration.
#' @param disease.prev prevalence of the binary outcome.
#' @param genotype a vector that represents the exposure data.
#' @param subject.effect.data subject effect data, reflects the heterogenity in baseline disease risk.
#' @param geno.OR odds ratio related to the 'at risk' genotype.
#' @return a binary vector that represents the phenotype data.
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # Generate data for a binary SNP with a MAF of 0.1
#' geno.elts <- sim.geno.data(num.obs=10000, geno.model=0, MAF=0.1)
#' geno.data <- geno.elts$genotype
#' 
#' # simulate subject effect data
#' s.effect <- sim.subject.data(num.obs=10000, baseline.OR=10)
#' 
#' # Generate binary phenotype statuses using the exposure data obtained above and an OR of 1.5 for the 'at risk' genotype
#' # and a prevalence of 0.1 for the binary outcome
#' pheno.data <- sim.pheno.bin.G(num.obs=10000,disease.prev=0.1, genotype=geno.data, subject.effect.data=s.effect, geno.OR=1.5)
#' 
#' }
#' 
sim.pheno.bin.G <- function(num.obs=10000, disease.prev=0.1, genotype=NULL, subject.effect.data=NULL, geno.OR=1.5){
  
  # IF GENOTYPE AND SUBJECT EFFECT DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
  if(is.null(genotype)){
    cat("\n\n ALERT!\n")
    cat(" No genotype data found.\n")
    cat(" Check the argument 'genotype'\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  if(is.null(subject.effect.data)){
    cat("\n\n ALERT!\n")
    cat(" No baseline effect data found.\n")
    cat(" Check the argument 'subject.effect.data'\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  numobs <- num.obs
  pheno.prev <- disease.prev
  genodata <- genotype
  s.efkt.data <- subject.effect.data
  geno.odds <- geno.OR
 
  # GET THE ALPHA AND BETA VALUES
  alpha <- log(pheno.prev/(1-pheno.prev))
  beta <- log(geno.odds)
  
  # GENERATE THE LINEAR PREDICTOR
  lp <- alpha + (beta*genodata) + s.efkt.data

  # GET 'mu' THE PROBABILITY OF DISEASE THROUGH LOGISTIC TRANSFORMATION
  mu <- exp(lp)/(1 + exp(lp))

  # GENERATE THE PHENOTYPE DATA AND RETURN IT AS A DATAFRAME
  phenotype <- rbinom(numobs,1,mu)
  return(phenotype)
}

