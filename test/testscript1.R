# this script test 18 scenarios, for a quantitative outcome, to compare the results across differents setting 
# and ensure there is no 'weird' results due to some bug not spotted by R checks.
#
# the 18 scenarios are divided into 6 classes of checks as follow (starting from scenario 1 to 18):
# 1st class: 3 minor allele frequencies; 0.15, 0.25 and 0.4
# 2nd class: 3 effect sizes related to the 'at risk' genetic exposure; 1.2, 2 and 3.5
# 3rd class: 2 genetic model; binary (0) and additive (1)
# 4th class: 3 sensitivity values for the outcome; 0.65, 8, 1
# 5th class: 3 specificity values for the outcome; 0.65, 8, 1
# 6th class: 4 start number of subjects ; 700, 800, 900, 1000
#
# apart from these parameters all the other parameters are set to the same values across all 18 scenarios

# load the library required for an analysis with a main effect environmental exposure
library(ESPRESSO.G)

# set the working directory to the 'test' directory
setwd(paste(getwd(), "/test", sep=""))

# load the input parameter tables
simulation.params <- read.table("simulation.params1.csv", header=T, sep=";")
pheno.params <- read.table("pheno.params1.csv", header=T, sep=";")
env.params <- read.table("geno.params1.csv", header=T, sep=";")

# call the function that runs a full ESPRESSO analysis and set it to run all 18 scenarios
run.espresso.G(simulation.params, pheno.params, env.params, scenarios2run=c(1:18))
