#' 
#' @title Fits a Generalized Linear Model
#' @description This function fits a conventional unconditional logistic regression model 
#' on the observed outcome and exposure data. In a Typical ESPRESSO the function is called
#' at each iteration and the estimates stored.
#' @param pheno.model outcome type: binary=0 and normal=1.
#' @param observed.data a dataframe that contains the observerd outcome and covariate data.
#' @return a vector that holds the beta, standard-error and z-statistic of each of the covariates.
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # load a table that contains observed quantitative outcome and binary SNP data
#' data(observed.data.G)
#'
#' # run a logistic regression (outcome is quantitative)
#' glm.estimates <- glm.analysis.G(pheno.model=1, observed.data.G)
#' 
#' }
#' 
glm.analysis.G <- function(pheno.model=0, observed.data=NULL){
  
  if(is.null(observed.data)){
			 cat("\n\n ALERT!\n")
			 cat(" No data found.\n")
			 cat(" Check the argument 'observed.data'\n")
			 stop(" End of process!\n\n", call.=FALSE)
		} 

  # BINARY OUTCOME
  if(pheno.model == 0){
	   # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
	   mod.glm <- glm(phenotype ~ 1 + genotype, family=binomial, data=observed.data)
	   mod.sum <- summary(mod.glm)
  }
  
  # QUANTITATIVE OUTCOME     
  if(pheno.model == 1){
	    # FIT A GLM FOR A GAUSSIAN OUTCOME
	    mod.glm <- glm(phenotype ~ 1 + genotype, family=gaussian, data=observed.data)
	    mod.sum <- summary(mod.glm)     
  }
  
	 beta.value <- mod.sum$coefficients[2,1]
	 se.value <- mod.sum$coefficients[2,2]
	 z.value <- mod.sum$coefficients[2,3]
	 
  # RETURN A VECTOR
  return(list(beta=beta.value, se=se.value, z=z.value))
}

