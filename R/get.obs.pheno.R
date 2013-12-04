#' 
#' @title Generates observed outcome data
#' @description Adds a set level of error to error free binary or quantitative data (the true phenotype data) 
#' to obtain data with a larger variance (the observed phenotype data).
#' @param phenotype outcome status.
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1 or uniform=2.
#' @param pheno.error misclassification rates: 1-sensitivity and 1-specificity
#' @param pheno.reliability reliability of the assessment for a quantitative outcome.
#' @return A dataframe containing:
#' \code{true.phenotype} the error free outcome data (true data).
#' \code{observed.phenotype} the true outcome data with some added error (observed data).
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # Example 1: a binary phenotype determined by a binary SNP
#' # Generate data for a binary SNP with a MAF of 0.1
#' geno.elts <- sim.geno.data(num.obs=10000, geno.model=0, MAF=0.1)
#' geno.data <- geno.elts$genotype
#' # generate subject effect data (i.e. heterogeneity in baseline risk) with a baseline OR of 10
#' s.effect <- sim.subject.data(num.obs=10000, baseline.OR=10)
#' # generate the 'true' (error free) phenotype data with OR of the exposure = 1.5
#' pheno.data <- sim.pheno.bin.G(num.obs=10000,disease.prev=0.1, genotype=geno.data, subject.effect.data=s.effect, geno.OR=1.5)
#' # get the observed phenotype data if the sensitivity and specificity of the assessment of the phenotype are both 0.95
#' obs.pheno.data <- get.obs.pheno(phenotype=pheno.data, pheno.model=0, pheno.error=c(0.05,0.05))
#' 
#' # Example 2: a quantitative phenotype determined by a binary SNP
#' # Generate data for a binary SNP with a MAF of 0.1
#' geno.elts <- sim.geno.data(num.obs=10000, geno.model=0, MAF=0.1)
#' geno.data <- geno.elts$genotype
#' # Generate quantitative phenotype statuses using the exposure data obtained above and an effect size of 0.25 
#' # for the 'at risk' genotype
#' pheno.data <- sim.pheno.qtl.G(num.subjects=10000, genotype=geno.data, geno.efkt=0.25)
#' # get the observed phenotype data if reliability of the assessment of the phenotype is 0.8
#' obs.pheno.data <- get.obs.pheno(phenotype=pheno.data, pheno.model=1, pheno.reliability=0.8)
#' 
#' }
#'
get.obs.pheno <- function (phenotype=NULL, pheno.model=0, pheno.error=c(0.05,0.05), pheno.reliability=0.9){ 
  
  if(is.null(phenotype)){
    cat("\n\n ALERT!\n")
    cat(" No phenotype data found.\n")
    cat(" Check the argument 'phenotype'\n")
    stop(" End of process!\n\n", call.=FALSE)
  } 
  if(is.null(pheno.model)){
    cat("\n\n ALERT!\n")
    cat(" No outcome  model provided\n")
    cat(" Check the argument 'pheno.model'\n")
    stop(" End of process!\n\n", call.=FALSE)
  } 
  
  true.phenotype <- phenotype
  misclass.rate.1.to.0 <- pheno.error[1]
  misclass.rate.0.to.1 <- pheno.error[2]	
  
  # GET THE OBSERVED OUTCOME DATA
  if(pheno.model==0){ # IF THE OUTCOME IS BINARY
    
    observed.phenotype <- misclassify(true.phenotype, misclass.rate.1.to.0, misclass.rate.0.to.1)
    
  }else{ # IF THE OUTCOME IS CONTINUOUS NORMAL
    
    # USE THE RELIABITLITY OF PHENOTYPE ASSESSMENT TO COMPUTE THE VARIANCE OF MEASURED PHENOTYPES.
    # RELIABITY = (VAR.true.data/VAR.obs.data) AND VAR.obs.data = (VAR.true.data + VAR.measurement)
    # IT FOLLOWS THAT VAR.true.data + VAR.of.estimate = VAR.true.data/RELIABILITY AND THEN:
    # VAR.measurement = (VAR.true.data/RELIABILITY) - VAR.true.data
    var.m <- (1^2/pheno.reliability)-(1^2) # standardized SD (SD = 1)
    
    # GENERATE THE NORMALLY DISTRIBUTED ERROR USING THE ABOVE COMPUTED VARIANCE
    num.obs <- length(phenotype)
    pheno.error <- rnorm(num.obs,0,sqrt(var.m))
    
    # ADD THE ERROR TO ORIGINAL PHENOTYPES TO GENERATE THE OBSERVED PHENOTYPE DATA
    observed.phenotype <- true.phenotype + pheno.error
  } 
  
  # RETURN THE TRUE AND OBSERVED PHENOTYPE DATA AS A DATAFRAME
  df <- data.frame(true.phenotype, observed.phenotype)
}