# pipelines

`pipelines` is a R package for performing different sorts of exploratory analyses on omics data sets.

## Installation

1. Install `devtools` from CRAN with `install.packages("devtools")`.

2. Install/Update `pipelines` with `devtools::install_github("erdto/TPDT")`

After that you load it as usual with `library(pipelines)`

## Demo

To give a short example, we will simulate some metabolite data, 3 groups and 2 covariates that should 
be controlled for.

```r
library(pipelines)
metdata <- matrix(rnorm(20 * 100, sd = 5), nrow = 100)
groups <- factor(sample(c("N", "FL", "DLBCL"), size = 100, replace = TRUE))
covariates <- data.frame(C_Sex = factor(sample(0:1, size = 100, replace = TRUE)),
                         C_Age_at_diagnosis = as.numeric(sample(20:60, 100, replace = TRUE)))

# run it
pipeline_pairwise(data = metdata, groups = groups, 
                  compname = "demo_pairwise",
                  vargroups = vargroups,
                  covars = covariates)

# The output is saved in an excel file and plots in pdfs in the working directory.
```

For further help, look into the vignette or the help pages with ?pipeline_pairwise.

## Volcano Plot
As a visual summary, among the saved pdfs are volcano plots for each combination of groups.
The observations beyond the cutoff lines in x and y direction can be considered worth to further investigate.

![plot1](figure/readme_volcano.pdf) 
