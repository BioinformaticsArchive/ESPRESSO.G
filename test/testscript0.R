# this script tests 18 scenarios, for a binary outcome, to compare the results across differents setting 
# and ensure there is no 'weird' results due to some bug not spotted by R checks.
#
# the 18 scenarios are divided into 6 classes of checks, as follows (starting from scenario 1 to 18):
# 1st class: 3 minor allele frequencies; 0.1, 0.2 and 0.3 (set in 'geno.param0s' file)
# 2nd class: 3 OR related to the 'at risk' genetic exposure; 1.2, 2 and 3.5 (set in 'geno.params0' file)
# 3rd class: 3 Baseline ORs (heterogeneity in baseline risk of disease); 8, 12 and 16 (set in 'simulation.params0' file)
# 4th class: 3 sensitivity values for the SNP alleles measurement; 0.65, 0.8, 0.95 (set in 'geno.params0' file)
# 5th class: 3 specificity values for the SNP alleles measurement; 0.65, 0.8, 0.95 (set in 'geno.params0' file)
# 6th class: 3 start number of cases/controls ; 1000/4000, 2000/8000 and 2000/4000 (set in 'simulation.params0' file)
#
# apart from these parameters all the other parameters are set to the same values across all 18 scenarios

# load the library required for an analysis with a main effect environmental exposure
library(ESPRESSO.G)

# set the working directory to the 'test' directory
#setwd(paste(getwd(), "/test", sep=""))

# load the input parameter tables
simulation.params <- read.table("simulation.params0.csv", header=T, sep=";")
pheno.params <- read.table("pheno.params0.csv", header=T, sep=";")
geno.params <- read.table("geno.params0.csv", header=T, sep=";")

# call the function that runs a full ESPRESSO analysis and set it to run all 18 scenarios
run.espresso.G(simulation.params, pheno.params, geno.params, scenarios2run=c(1:18))
