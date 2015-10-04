#library(xlsx)
#install.packages("ggplot2")
#install.packages("XLConnect")
options(java.parameters = "-Xmx1024m")
# library(ggplot2)
# library(XLConnect)
# library(checkmate)
# library(parallel)

# load utility functions
# source(file.path(getwd(), "basic_pipelines", "utilities_pipeline.R"))



#' Pairwise group comparisons of metabolomics data.
#'
#' Computes pairwise t-tests (or linear regressions when covariates present) or fisher-tests, of the difference of all groups in all metabolites given in "data".
#' 
#' @param data  [data.frame] with metabolites as columns.
#' @param groups  [vector] grouping variable with length matching the number of rows of \code{data}
#' @param covars [data.frame] with variables to be controlled for.
#' @param compname  [character] prefix for name of output files (excel file, boxplots, volcano plots).
#' @param vargroups [data.frame] with columns being grouping variables of the metabolites in \code{data}. 
#' Used for computing aggregated z-scores for the groups to be compared instead of comparing the single metabolites.
#' @param test [character] what test is to be performed? Linear regression (if no covariates provided, equivalent to t-test) or
#' fishers exact test for comparing the number of NA (below threshold concentration) vs. non-NA values per metabolite.
#' @param ratios [logical] Should ratios of the metabolites be computed an tested? Will quickly lead to long runtimes.
#' @param ratioop [character] Should ratios be computed by "subtract"ing the metabolite values (for logarithmic data) or by "divid"ing?
#' @param p_cutoff [numeric] Only metabolites with below cut-off p-value will be shown in boxplots. Default is no cut-off.
#' @param repetitions [numeric] Number of repetitions for computing resampling-corrected p-values. Usually, at least 200 are recommended,
#' though higher accuracy can only obtain with greater numbers of repetitions.
#' @param boxplot_args deprecated.
#' @param volc_args deprecated.
#' @param no_plots [logical] Should all plots be omitted? This will speed up the function considerably.
#' @param verbose [logical] Should progress be printed?
#' @param ncores [numeric] number of CPUs to be used. If >1, computations will be parallelized.
#' @return Returns NULL, all output (excel file, pdf with boxplots and volcano plots) is saved into the working directory.
#' @export

pipeline_continuous <- function(data, predictors, compname, 
                                vargroups=data.frame(),
                                cor_method="pearson",
                                test = c("linreg", "fisher"),
                                covars=data.frame(), 
                                ratios=FALSE, 
                                ratioop="subtract", 
                                ratio_cutoff=.1, 
                                repetitions=0,
                                volc_args =  set_volc_args(),
                                no_plots = FALSE,
                                ncores = 1) {
  
  test <- match.arg(test)
  checkDataFrame(data, types = "numeric", all.missing = TRUE)
  vars <- data
  varnames <-  names(data)
  covarnames <- names(covars)
  prednames <- names(predictors)
  predictors <- as.data.frame(predictors)
  
  vargroups <- na.exclude(vargroups)
  vncols <- ncol(vargroups)
  
  if(ratios == TRUE & no_plots == FALSE)
    warning(paste("Warning, ratios = TRUE and no_plots = FALSE, this could take a while!", "\n",
                  "Using p_cutoff argument is recommended when ratios is on."))
  
  # We use 'data' inside the loops later on, so we save the original in
  # 'data1' to call 'use_agg_z_df' with at the beginning of the vargroup-loops
  
  #   data1 <- data
  # changed, we use vars in the loop and use data to get agg_z scores at the beginning of the
  # vargroup iterations.
  
  if(ratios==TRUE){
    data <- add_ratios_to_df(df = data, names = varnames)
  }
  
  
  #   results_list <- vector("list", vncols + 1)
  
  #create file to save workbook to
  file <- paste(getwd(), "/", compname, ".xls", sep="")
  
  # alt: 
  #wb <- createWorkbook(type="xls")
  wb <- XLConnect::loadWorkbook(file, create=TRUE)
  
  for(vloop in 0:vncols) {
    
    if(vloop != 0) {
      cat("\r", "vargroup-loop:", vloop, "out of", vncols) #, "\n") 
      
      # y != 0 means we're in a vargroup-loop, new data matrix with agg z scores will be used
      vars <- use_agg_z_df(df = data, vargrp = vargroups[ , vloop])
      
      varcomb <- unique(vargroups[ , vloop])
      #       newdata <- data.frame(matrix(NA, nrow=nrow(data), ncol=length(varcomb)))
      #       names(newdata) <- varcomb
      
      if(ratios == TRUE) {
        varnames <- names(data)
        data <- add_ratios_to_df(df = data, names = varnames)
      }
      
    }
    
    results <- as.data.frame(matrix(NA, nrow=NCOL(vars),ncol= 10))
    names(results) <- c("inputvar","predictor", "N", "beta", "pvalue", "p_fdr", "p_bonf", "p_sampling", "CIlower","CIupper") 
    results$inputvar <- colnames(vars)
    results_list <- vector("list", length(prednames))
    
    for(i in 1:length(prednames)) {
      
      # open scatterplot pdf
      #       pdf(file=paste(compname, "_",prednames[i], "_", "scatterplots", 
      #                      if(vloop!=0){paste("_", colnames(vargroups)[vloop])}, ".pdf", sep=""), onefile=T)
      
      
      #initialize matrix results where results will be saved and which is written into excel sheets later
      results <- as.data.frame(matrix(NA, nrow=NCOL(vars),ncol= 10))
      names(results) <- c("inputvar", "predictor", "N", 
                          "beta", "pvalue", "p_fdr", "p_bonf", "p_sampling", "CIlower","CIupper") 
      results$inputvar <- colnames(vars)
      results$predictor <- prednames[i]
      
      
      
      number_cases <- sapply(vars, function(x) sum(!is.na(x)))
      results$N <- number_cases
      indeces <- (number_cases >= 3)
      notnas <- names(number_cases)[indeces]
      # skip vars where all rows are NA
      vars <- vars[ , colnames(vars) %in% notnas]
      
      
      #         predictors[[prednames[i]]]
      
      
      output <-
        lapply(names(vars), function(j) {
          get_row2(X = data.frame(y = vars[,j], x = predictors[,i]),
                   name = prednames[i], covars, test = test)
         })
      
      results <- data.frame(inputvar = names(vars), dplyr::rbind_all(output))
#         do.call(rbind, output)
      
      if(repetitions > 0) {
        cat("\n", "resampling...", "\n")
        
        # resampling: vector for correction
        p_rand <- numeric(repetitions)
        
        # do resampling in chunks of parallel threads
        chunks <- suppressWarnings( split(1:repetitions, 1:ncores) )
        
        
        if(.Platform$OS.type != "windows" || ncores == 1) {
          
          p_rand <- unlist(mclapply(chunks, function(repchunk) {
            
            sapply(repchunk, function(rep) {
              min(sapply(vars, function(var) 
                get_row2(X = data.frame(y = sample(var),
                                    x = predictors[, i]),
                     covars = covars, test = test, name = prednames[i])$p),
                na.rm = TRUE)
            })
        }, mc.cores = ncores))
        
        } else {
          
          cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
          on.exit(parallel::stopCluster(cl))
          
          # functions called by compare_grps have to be exported explicitly
          parallel::clusterExport(cl = cl,  varlist = c("get_row2"), envir = environment())
          
          # if your function does not call "require(ggplot2)"
          #           parallel::clusterEvalQ(cl, library(ggplot2))
          
          p_rand <- unlist(parallel::parLapply(cl = cl, X = chunks, function(X) {
            sapply(X, function(rep) {
              min(sapply(vars, function(var) 
                get_row2(X = data.frame(y = sample(var),
                                        x = predictors[, i]),
                         covars = covars, test = test, name = prednames[i])$p),
                na.rm = TRUE)
            })
          }))
        }
        
        # add correted pvals to results data.frame
        results$p_sampling <- sapply(results[,"pvalue"], function(x) sum(unlist(p_rand) < x) / repetitions)
      } # here resampling is done over    
     
      # after resampling,
      # add sampling p-values to dataframe
      results$p_fdr <- p.adjust(results$pvalue, method="fdr")
      results$p_bonf <- p.adjust(results$pvalue, method="bonferroni")
      results_list[[i]] <- results
      
      # vloop=0 -> metabolite loop
      
      if(ratios && vloop==0) {
       results <- compute_pgain(R = results)
      }
      # scatterplots
      if(all(c(no_plots == FALSE, test=="linreg", 
               length(results$beta) > 1, sum(!is.na(results$beta)) > 2))) {
        pdf(paste(compname, "_", prednames[i], "_", "scatterplot", 
                  if(vloop!=0){paste("_", colnames(vargroups)[vloop])},".pdf", sep=""))
        
        sapply(names(vars), function(j) print_scatterplot(y = vars[, j], 
                                                         predictor = predictors[, i], 
                                                         covariates = covariates, 
                                                         name = j, predname = prednames[i]))  
        dev.off()
      }
               
      #volc
      if(all(c(no_plots == FALSE, test=="linreg", length(results$beta) > 1, 
               sum(!is.na(results$beta)) > 2))) {
        
        pdf(paste(compname, "_", prednames[i], "_", "volcano", 
                  if(vloop!=0){paste("_", colnames(vargroups)[vloop])},".pdf", sep=""))
        volc_args
        print_volc(datframe = results, volc_args = volc_args)
        dev.off()
      } else {
        paste("Omitted volcano for predictor: ", i)
      }
      
    } # here loop over predictors ends
    
    # write results to excel file
    RR <- do.call(rbind, results_list)
    SN <- paste(compname, if(vloop!=0){paste("_",colnames(vargroups)[vloop])}, sep="")
    sheet <- XLConnect::createSheet(wb, name=SN)
    
    # save data sorted by pvalue
    XLConnect::writeWorksheet(wb, RR[with(RR, order(pvalue)), ], sheet = SN)
    # before: just RR
    
  } # here loop over vargroups ends
  XLConnect::saveWorkbook(wb, file=file)
  cat("\n", "done.")
}