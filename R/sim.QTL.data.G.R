#' 
#' @title Simulates subjects for a continuous outcome
#' @description Generates the specified number of subjects for a quantitative outcome.
#' @param numsubjects number of subjects to simulate.
#' @param ph.mean statistical mean
#' @param ph.sd standard deviation
#' @param geno.model genetic model; binary=0 and additive=1.
#' @param MAF minor allele frequencies of the genetic variant.
#' @param geno.efkt effect size of the 'at risk' genotype.
#' @param pheno.reliability reliability of the assessment for a quantitative outcome.
#' @return a matrix that holds the outcome (\code{phenotype}) and exposure (\code{genotype}) data 
#' and the SNP alleles used to construct the genotypes.
#' @keywords internal
#' @author Gaye A.
#' 
sim.QTL.data.G <- function(numsubjects=10000, ph.mean=0, ph.sd=1, geno.model=0, MAF=0.1, geno.efkt=0.25, pheno.reliability=0.9){
  
	 geno.mod <- geno.model
   geno.maf <- MAF
   geno.efsize <- geno.efkt
   pheno.rel <- pheno.reliability
   
   # GENERATE THE TRUE GENOTYPE DATA
	 geno.data <- sim.geno.data(num.obs=numsubjects, geno.model=geno.mod, MAF=geno.maf)
	 allele.A <- geno.data$allele.A
	 allele.B <- geno.data$allele.B
	 genotype <- geno.data$genotype

   # GENERATE THRUE OUTCOME DATA
   genodata <- genotype
   pheno.data <- sim.pheno.qtl.G(num.subjects=numsubjects, pheno.mean=ph.mean, pheno.sd=ph.sd, genotype=genodata, geno.efkt=geno.efsize)
   true.phenotype <- pheno.data
   
	 # GENERATE THE OBSERVED OUTCOME DATA 
	 obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=1, pheno.sd=ph.sd, pheno.reliability=pheno.rel)
	 phenotype <- obs.phenotype$observed.phenotype
   
   # STORE THE GENERATED OBSERVED OUTCOME AND TRUE EXPOSURE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(phenotype,genotype,allele.A,allele.B)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id", "phenotype", "genotype", "allele.A", "allele.B")
   mm <- data.frame(sim.matrix)
   
}

