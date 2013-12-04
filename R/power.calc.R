#' 
#' @title Calculates the empirical and theoretical power
#' @description The function determines the empirical and theoretical power. 
#' The empirical power is the proportion of simulations in which 
#' the z-statistic for the parameter of interest exceeds the z-statistic 
#' for the desured level if statistical significance. 
#' The theoretical power is the power of the study.
#' @param pval cut-off p-value defining statistical significance.
#' @param z.values z-statistic of the determinant.
#' @param mean.model.z mean z-statistic of the environmental determinant.
#' @return a list that contains the computed empirical power and theoretical power.
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # generate the arguments to pass on to the function (typically the output of the function glm.analysis.G)
#' # load a table of quantitative outcome and binary SNP data
#' data(observed.data.G)
#'
#' # run a logistic regression (outcome is binary)
#' glm.estimates <- glm.analysis.G(pheno.model=1, observed.data.G)
#' beta <- glm.estimates$beta
#' se <- glm.estimates$se
#' modz <- beta/se
#' zvals <-  glm.estimates$z
#' 
#' # use the function; in a typical ESPRESSO process 'modz' would be the ratio of the average beta over the average
#' # se across the stored beta and se values (after each run the estimates of the glm analysis are stored) and 'zvals'
#' # would be the vector of z values obtained after the iterations.
#' power.values <- power.calc(z.values=zvals, mean.model.z=modz)
#'}
#'
power.calc <- function(pval=1e-04, z.values=NULL, mean.model.z=NULL){
  
  if(is.null(z.values)){
    cat("\n\n ALERT!\n")
    cat(" No z-statistics found\n")
    cat(" Check the argument 'z.values'.\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  if(is.null(mean.model.z)){
    cat("\n\n ALERT!\n")
    cat(" The argument 'mean.model.z' is set to NULL.\n")
    cat(" This argument should be the ratio 'mean.beta/mean.se'.\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  # CALCULATE Z STATISTIC THRESHOLD FOR DESIRED P-VALUE 
  z.pval <- qnorm(1-pval/2)
  
  # GET EMPIRICAL POWER: THE PROPORTION OF SIMULATIONS IN WHICH THE 
  # Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC 
  # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
  empirical.power <- round(mean((z.values > z.pval), na.rm=TRUE),3)
  
  # GET THE MODELLED POWER
  modelled.power <- pnorm(mean.model.z-z.pval)
  
  return(list(empirical=empirical.power, modelled=modelled.power))
}
