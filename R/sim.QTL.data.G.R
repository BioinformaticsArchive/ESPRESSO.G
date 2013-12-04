#' 
#' @title Simulates subjects for a continuous outcome
#' @description Generates the specified number of subjects for a quantitative outcome.
#' @param numsubjects number of subjects to simulate.
#' @param geno.model genetic model; binary=0 and additive=1.
#' @param MAF minor allele frequencies of the genetic variant.
#' @param geno.efkt effect size of the 'at risk' genotype.
#' @return a matrix that holds the outcome (\code{phenotype}) and exposure (\code{genotype}) data 
#' and the SNP alleles used to construct the genotypes.
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # Example 1: generate 5000 subjects for a binary genetic variant with a MAF of 0.1 and an effect size of 0.25
#' sim.matrix <- sim.QTL.data.G(numsubjects=5000, geno.model=0, MAF=0.1, geno.efkt=0.25)
#'
#' # Example 2: generate 5000 subjects for an additive genetic variant with a  MAF of 0.1 and an effect size of 0.25
#' sim.matrix <- sim.QTL.data.G(numsubjects=5000, geno.model=1, MAF=0.1, geno.efkt=0.25)
#' 
#' }
#' 
sim.QTL.data.G <- function(numsubjects=10000, geno.model=0, MAF=0.1, geno.efkt=0.25){
  
   num.obs <- numsubjects
	   
   # GENERATE THE TRUE GENOTYPE DATA
			geno.data <- sim.geno.data(num.obs, geno.model, MAF)
			allele.A <- geno.data$allele.A
			allele.B <- geno.data$allele.B
			genotype <- geno.data$genotype

   # GENERATE OUTCOME DATA
   pheno.data <- sim.pheno.qtl.G(num.obs, genotype, geno.efkt)
   phenotype <- pheno.data

   # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(phenotype,genotype,allele.A,allele.B)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id", "phenotype", "genotype", "allele.A", "allele.B")
   mm <- data.frame(sim.matrix)
}

