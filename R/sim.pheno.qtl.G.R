#' 
#' @title Simulates continuous outcome data
#' @description The function uses the effects data of the genetic determinant to construct a linear predictor(LP). 
#' The outcome is a normally distributed variable generated with a mean equal to the LP and a standard 
#' deviation of 1.
#' @param num.subjects number of subjects to generate.
#' @param pheno.mean statistical mean
#' @param pheno.sd standard deviation
#' @param genotype a vector that represents the exposure data
#' @param geno.efkt effect size of related to the 'at risk' allele.
#' @return a binary vector that represents the phenotype data.
#' @keywords internal
#' @author Gaye A.
#'  
sim.pheno.qtl.G <- function(num.subjects=10000, pheno.mean=0, pheno.sd=1, genotype=NULL, geno.efkt=0.25){
  
   # IF GENOTYPE DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
   if(is.null(genotype)){
			  cat("\n\n ALERT!\n")
			  cat(" No genotype data found.\n")
			  cat(" Check the argument 'genotype'\n")
			  stop(" End of process!\n\n", call.=FALSE)
		}
   
   numobs <- num.subjects
   genodata <- genotype
   geno.efsize <- geno.efkt

   # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   alpha <- pheno.mean
   beta <-	geno.efkt
   num.obs <- num.subjects

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (beta*genotype)

   # GENERATE THE TRUE PHENOTYPE DATA TO RETURN
   phenotype <- rnorm(num.obs,lp,pheno.sd)
   return(phenotype)
}

