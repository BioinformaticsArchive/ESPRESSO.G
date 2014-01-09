#' 
#' @title Simulates case and controls
#' @description Generates affected and non-affected subjects until the set sample size is achieved
#' @param block.size number of observations to generate per iteration until the specified number 
#' of cases and controls is achieved.
#' @param numcases number of cases to generate.
#' @param numcontrols number of controls to generate.
#' @param allowed.sample.size maximum number of observations allowed i.e. the total size of the 
#' population to sample from.
#' @param disease.prev prevalence of the binary outcome.
#' @param geno.model genetic model; binary=0 and and additive=1.
#' @param MAF minor allele frequency of the genetic variant.
#' @param geno.OR odds ratio of the genetic determinant.
#' @param baseline.OR baseline odds ratio for subject on 95 percent population centile versus 5 percentile. 
#' This parameter reflects the heterogeneity in disease risk arising from determinates that have not been 
#' measured or have not been included in the model.
#' @return a list which holds a matrix, \code{data}, that contains the phenotype and genotype statuses
#' and an integer, \code{allowed.sample.size.exceeded}, which tells if the maximum population size has been 
#' exceeded.
#' @export
#' @author Amadou Gaye
#' @examples {
#' 
#' # number of cases 
#' cases <- 2000
#' # number of controls
#' controls <- 8000
#'
#' # Example 1: generate cases and controls untill the set number of cases and controls (sample size) is achieved for
#' # a binary SNP with a MAF of 0.1 and an OR of 1.5;
#' # the heterogeneity in baseline risk of disease (baseline.OR) is 10.
#' sim.matrix <- sim.CC.data.G(block.size=20000,numcases=cases,numcontrols=controls,allowed.sample.size=20000000,
#'                           disease.prev=0.1, geno.model=0, MAF=0.1, geno.OR=1.5, baseline.OR=10)
#' data.generated <- sim.matrix$data
#'
#' # Example 2: generate cases and controls untill the set number of cases and controls (sample size) is achieved for
#' # an additive SNP with a MAF of 0.1 and an OR of 1.5;
#' # the heterogeneity in baseline risk of disease (baseline.OR) is 10.
#' sim.matrix <- sim.CC.data.G(block.size=20000,numcases=cases,numcontrols=controls,allowed.sample.size=20000000,
#'                           disease.prev=0.1, geno.model=1, MAF=0.1, geno.OR=1.5, baseline.OR=10)
#' data.generated <- sim.matrix$data
#' 
#' }
#' 
sim.CC.data.G <- function(block.size=20000, numcases=2000, numcontrols=8000, allowed.sample.size=20000000, 
                        disease.prev=0.1, geno.model=0, MAF=0.1, geno.OR=1.5, baseline.OR=12.36){
  
   numobs <- block.size
   ncases <- numcases
   ncontrols <- numcontrols
   max.pop.size <- allowed.sample.size
   pheno.prev <- disease.prev
   geno.mod <- geno.model
   geno.maf <- MAF
   geno.odds <- geno.OR
   baseline.odds <- baseline.OR
  
   # SET UP ZEROED COUNT VECTORS TO DETERMINE WHEN ENOUGH CASES AND CONTROLS HAVE BEEN GENERATED
   complete <- 0
   complete.absolute <- 0
   cases.complete <- 0
   controls.complete <- 0
   block <- 0

   # SET UP A MATRIX TO STORE THE GENERATED DATA
   sim.matrix <- matrix(numeric(0), ncol=4)

   # SET LOOP COUNTER
   numloops <- 0

   # LOOP UNTIL THE SET NUMBER OF CASES AND OR CONTROLS IS ACHIEVED OR THE 
   # THE SET POPULATION SIZE TO SAMPLE FROM IS REACHED
   while(complete==0 && complete.absolute==0){

     # GENERATE THE TRUE GENOTYPE DATA
     geno.data <- sim.geno.data(num.obs=numobs, geno.model=geno.mod, MAF=geno.maf)
     allele.A <- geno.data$allele.A
     allele.B <- geno.data$allele.B
     genotype <- geno.data$genotype
     
     # GENERATE SUBJECT EFFECT DATA THAT REFLECTS BASELINE RISK: 
     # NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR WITH APPROPRIATE 
     # VARIANCE ON SCALE OF LOG-ODDS
     subject.effect.data <- sim.subject.data(num.obs=numobs, baseline.OR=baseline.odds)

     # GENERATE THE TRUE OUTCOME DATA
     genodata <- genotype
     s.efkt.data <- subject.effect.data
     pheno.data <- sim.pheno.bin.G(num.obs=numobs, disease.prev=pheno.prev, genotype=geno.data$genotype, 
                                   subject.effect.data=s.efkt.data, geno.OR=geno.odds)
     phenotype <- pheno.data
     
     # STORE THE TRUE OUTCOME, GENETIC AND ALLELE DATA IN AN OUTPUT MATRIX 
     # WHERE EACH ROW HOLDS THE RECORDS OF ONE INDIVUDAL
     sim.matrix.temp <- cbind(phenotype,genotype,allele.A,allele.B)
     
     # UPDATE THE MATRIX THAT HOLDS ALL THE DATA GENERATED SO FAR, AFTER EACH LOOP
     sim.matrix <- rbind(sim.matrix, sim.matrix.temp)
     
     # SELECT OUT CASES
     sim.matrix.cases <- sim.matrix[phenotype==1,]
     
     # SELECT OUT CONTROLS
     sim.matrix.controls <- sim.matrix[phenotype==0,]
     
     # COUNT THE NUMBER OF CASES AND CONTROLS THAT HAS BEEN GENERATED
     cases.simulated <- dim(sim.matrix.cases)[1]
     controls.simulated <- dim(sim.matrix.controls)[1]
     
     # TEST IF THERE ARE AT LEAST ENOUGH CASES ALREADY SIMULATED
     # IF THERE ARE, DEFINE THE CASE ELEMENT OF THE DATA MATRIX
     if(cases.simulated >= numcases)
     {
       sim.matrix.cases <- sim.matrix.cases[1:numcases,]
       cases.complete <- 1
     }
     
     # TEST IF THERE ARE AT LEAST ENOUGH CONTROLS ALREADY SIMULATED
     # IF THERE ARE, DEFINE THE CONTROL ELEMENT OF THE DATA MATRIX
     if(controls.simulated>=numcontrols)
     {
       sim.matrix.controls <- sim.matrix.controls[1:numcontrols,]
       controls.complete <- 1
     }
     
     # HAVE WE NOW GENERATED THE SET NUMBER OF CASES AND CONTROLS?
     complete <- cases.complete*controls.complete		
     
     # HAVE WE EXCEEDED THE TOTAL SAMPLE SIZE ALLOWED?
     complete.absolute <- (((block+1)*block.size)>=allowed.sample.size)
     if(complete.absolute==1) {sample.size.excess <- 1}else{sample.size.excess <- 0}
     
     # INCREMENT LOOP COUNTER
     numloops <- numloops + 1
   }

   # STACK FINAL DATA MATRIX WITH CASES FIRST
   sim.matrix <- rbind(sim.matrix.cases,sim.matrix.controls)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # NAME THE COLUMNS OF THE MATRIX AND RETURN IT AS A DATAFRAMEDATAFRAME
   colnames(sim.matrix) <- c("id", "phenotype", "genotype", "allele.A", "allele.B")
   mm <- list(data=data.frame(sim.matrix), allowed.sample.size.exceeded=sample.size.excess)
}

