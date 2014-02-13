#' 
#' @title Generates outcome and exposure data with some added error
#' @description Adds a specified level of error to exposure data to generate the 'observed' data.
#' See \code{get.obs.pheno} and \code{get.obs.geno} for details on how the observed outcome and exposure data are
#' obtained.
#' @param true.data input table of simulated data considered as error free.
#' @param geno.model Genetic model; 0 for binary and 1 for continuous
#' @param MAF minor allele frequency of the SNP (in ESPRESSO this is the frequency of the 'at risk' allele)
#' @param geno.error misclassification rates in the assessment of the SNP alleles: 1-sensitivity and 1-specificity
#' @return A matrix which contains the observed outcome and exposure data
#' @export
#' @author Gaye,A.
#' @examples {
#' 
#' # load a table of binary outcome and binary exposure (SNP) data
#' data(true.data.G)
#' 
#' # generate the 'observed' data by adding some error determined by the respective sensitivity and specificity 
#' # levels of the assessment of the exposure.
#' observed.data <- get.observed.data.G(true.data=true.data.G, geno.model=0, MAF=0.1, geno.error=c(0.05,0.05))
#' }
#' 
get.observed.data.G <- function(true.data=NULL, geno.model=0, MAF=0.1, geno.error=c(0.05,0.05)){               
  
    if(is.null(true.data)){
			   cat("\n\n ALERT!\n")
			   cat(" No data found.\n")
			   cat(" Check the argument 'true.data'\n")
			   stop(" End of process!\n\n", call.=FALSE)
		 	}

    sim.df <- true.data  
    geno.mod <- geno.model
    geno.maf <- MAF
    geno.err <- geno.error
    
    # GET THE OBSERVED GENOTYPES
	   true.genotype <- sim.df$genotype
	   obs.genotype <- get.obs.geno(allele.A=sim.df$allele.A, allele.B=sim.df$allele.B, 
                                  geno.model=geno.mod, MAF=geno.maf, geno.error=geno.err)
	   
	  # REPLACE THE TRUE GENOTYPE DATA BY THE NOW GENERATED OBSERVED GENOTYPES
	  # IN THE INITIAL MATRIX THAT HELD THE TRUE DATA
	  sim.df$genotype <- obs.genotype$observed.genotype
	  sim.df$allele.A <- obs.genotype$observed.allele.A
	  sim.df$allele.B <- obs.genotype$observed.allele.B
	   
    # RETURN THE MATRIX WHICH NOW CONTAINS ONLY THE OBSERVED DATA TO ANALYSE BY GLM
    colnames(sim.df) <- c("id", "phenotype", "genotype", "allele.A", "allele.B")
    return(sim.df)
    
}

