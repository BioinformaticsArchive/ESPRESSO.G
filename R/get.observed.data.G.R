#' 
#' @title Generates outcome and exposure data with some added error
#' @description Adds a specified level of error to outcome and exposure data to generate the 'observed' data.
#' See \code{get.obs.pheno} and \code{get.obs.geno} for details on how the observed outcome and exposure data are
#' obtained.
#' @param true.data input table of simulated data considered as error free.
#' @param geno.error misclassification rates in the assessment of the SNP alleles: 1-sensitivity and 1-specificity
#' @param geno.model Genetic model; 0 for binary and 1 for continuous
#' @param MAF minor allele frequency of the SNP (in ESPRESSO this is the frequency of the 'at risk' allele)
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1 or uniform=2.
#' @param pheno.error misclassification rates: 1-sensitivity and 1-specificity
#' @param pheno.reliability reliability of the assessment of a quantitative phenotype
#' @return A matrix which contains the observed outcome and exposure data
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # load a table of binary outcome and binary exposure (SNP) data
#' data(true.data.G)
#' 
#' # generate the 'observed' data by adding some error determined by the respective sensitivity and specificity 
#' # levels of the assessment of the outcome and the exposure.
#' observed.data <- get.observed.data.G(true.data=true.data.G, pheno.model=0, pheno.error=c(0.1,0.1), geno.model=0, 
#'                                      MAF=0.1, geno.error=c(0.05,0.05))
#' }
#' 
get.observed.data.G <- function(true.data=NULL, pheno.model=0, pheno.error=c(0.1,0.1), pheno.reliability=0.9, 
                                geno.model=0, MAF=0.1, geno.error=c(0.05,0.05)){               
  
    if(is.null(true.data)){
			   cat("\n\n ALERT!\n")
			   cat(" No data found.\n")
			   cat(" Check the argument 'true.data'\n")
			   stop(" End of process!\n\n", call.=FALSE)
		 	}
		 	
		 	sim.df <- true.data      

    # GET THE OBSERVED GENOTYPES
	   true.genotype <- sim.df$genotype
	   obs.genotype <- get.obs.geno(sim.df$allele.A,sim.df$allele.B, geno.model, MAF, geno.error)
	   
	   # REPLACE THE TRUE GENOTYPE DATA BY THE NOW GENERATED OBSERVED GENOTYPES
	   # IN THE INITIAL MATRIX THAT HELD THE TRUE DATA
	   sim.df$genotype <- obs.genotype$observed.genotype
	   sim.df$allele.A <- obs.genotype$observed.allele.A
	   sim.df$allele.B <- obs.genotype$observed.allele.B
	   
    # GET THE OBSERVED OUTCOME DATA
    true.phenotype <- sim.df$phenotype
    if(pheno.model==0){
      obs.phenotype <- get.obs.pheno(true.phenotype, pheno.model, pheno.error)
    }else{
      obs.phenotype <- get.obs.pheno(true.phenotype, pheno.model, pheno.reliability=pheno.reliability)
    }
    
    # REPLACE THE TRUE PHENOTYPE DATA BY THE NOW GENERATED OBSERVED PHENOTYPES
    sim.df$phenotype <- obs.phenotype$observed.phenotype
    
    # RETURN THE MATRIX WHICH NOW CONTAINS ONLY THE OBSERVED DATA TO ANALYSE BY GLM
    colnames(sim.df) <- c("id", "phenotype", "genotype", "allele.A", "allele.B")
    return(sim.df)
}

