# Demo

library(pipelines)

# Simulate data
metdata <- data.frame(matrix(rnorm(20 * 100, sd = 5), nrow = 100))
names(metdata) <- c("tryptophanBetaine", "glutamine", "tryptophan", "histidine", 
                 "leucine", "cholesterol", "phenylalanine", "creatinine", "lactate", 
                 "X4.hydroxyphenylacetate", "X3.hydroxybutyrate..BHBA.", "adenosine", 
                 "arabinose", "fructose", "mannose", "pyruvate", "uridine", "linoleate..18.2n6.", 
                 "allantoin", "arachidonate..20.4n6.")

# the grouping variable: All combinations of groups are tested pairwise
groups <- factor(sample(c("N", "FL", "DLBCL"), size = 100, replace = TRUE))

# covariables that are to be controlled for
covariates <- data.frame(C_Sex = factor(sample(0:1, size = 100, replace = TRUE)),
                         C_Age_at_diagnosis = as.numeric(sample(20:60, 100, replace = TRUE)))

# load variable groupings (Metabolites that can be grouped, for example into
# "Amino Acids" or "Glutamate Metabolism").
vargroups <- data.frame(Subpathway = sample(c("Xenobiotics", "Amino Acid", 
                      "Lipid", "Carbohydrate", "Nucleotide", 
                      "Energy", "Cofactors and Vitamins", 
                      "Peptide", NA_character_), 
                    size = 100, replace = TRUE))

##### Function call #########

# The function prints an excel file and several pdfs with the prefix
# "Example_pairwise" into the working directory

pipeline_pairwise(data = metdata, groups = groups,
                  compname = "Example_pairwise",
                  vargroups = vargroups, no_plots = TRUE,
                  covars = covariates)

# again with fisher's exact test,
# _fisher suffix will be added to the file name automatically.
pipeline_pairwise(data = metdata, groups = groups,
                  compname = "Example_pairwise",
                  vargroups = vargroups, no_plots = TRUE,
                  covars = covariates, test = "fisher")

##### pipeline_overall #########

# On the same data, the pipeline_overall can be applied. This compares
# all groups simultaneously with a F-test (can also control for covariates with
# regression). The F-test shows significant if at least one group differs from
# one of the rest.

pipeline_overall(data = metdata, groups = groups,
                 compname = "Example_overall",
                 vargroups = vargroups)

# now a fisher test
pipeline_overall(data = metdata, groups = groups,
                 refcat = "N", test = "fisher",
                 compname = "Example_overall",
                 vargroups = vargroups)


##### pipeline_continuous #########

# This function does not compare groups but regresses an continuous predictor on the metabolite
# data.
# For each metabolite, we will get a p-value for the effect of all variables in "predictors".
# In this case, we only use "age at diagnosis" as a predictor.
predictors <- data.frame(C_Age = realdata[, "C_Age_at_diagnosis__Y_"])

# You could also take some of the metabolites as predictors:
predictors <- data.frame(C_Age = covariates$C_Age_at_diagnosis,
                         glutamine = metdata$glutamine)

# We still can control i.e. for Sex
covars <- data.frame(C_Sex = covariates$C_Sex)

pipeline_continuous(data = metdata, 
                    predictors = predictors,
                    vargroups = vargroups,
                    compname = "demo_cont", 
                    covars = covars, 
                    repetitions = 0)

# Look into your working directory ( /R ) to find the volcano plots 
# and an excel file with the results.