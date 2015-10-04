### Utility function for pipeline_pairwise
compare_grps <- function(df, groups, vloop, compname, vargroups = data.frame(), covariates = data.frame(),
                         repetitions = 0, test = c("linreg", "fisher"), no_plots = FALSE,
                         boxplot_args, volc_args, verbose = TRUE) {
  compare <- levels(groups)
  # sanity check
  stopifnot(length(compare) == 2)
  
  if(no_plots) {
    boxplot_args <- set_boxplot_args()
    volc_args <- set_volc_args()
  }
 
  test <- match.arg(test)
  current_grp <- groups
  current_grp <- droplevels(current_grp)
  current_grp[which(!(groups %in% levels(groups)))] <- NA
  
  # resampling: vector for correction
  p_rand <- numeric(repetitions)
  
  # the following is done once:
  # clear rows of non-current-grp
  vars <- df
  vars[is.na(current_grp), ] <- NA
  
  # get sample sizes and indices of metabolites, that are not all na
  # fill in sample sizes
  index1 <- which(current_grp == compare[1])
  index2 <- which(current_grp == compare[2])
  number_cases <- sapply(vars, function(x) c(sum(!is.na(x[index1])),
                                             sum(!is.na(x[index2])) ) )
  notnas <- colnames(number_cases)[colSums(number_cases >= 2) == 2]
  
  #   if(test == "linreg") {
  #     # skip vars where all rows are NA
  #     vars <- vars[ , which(colnames(vars) %in% notnas)]
  #   }
  
  # now for reps...
  for(rloop in 0:repetitions) {
    
    # this once
    if(rloop==0) {
      
      # fill in results for output in excel sheet
      R <- data.frame(inputvar = names(df),
                      dplyr::rbind_all(
                        lapply(vars, function(v)
                          get_row1(y = v, group = current_grp, 
                                   covariates = covariates, test = test)))                            
      )
    }
    
    # now for rloop in 1:reps, we're resampling
    # so permutate grps and then get results
    if(rloop > 0) {
      if(rloop==1 && verbose) {
        cat("\n", "resampling...", "\n")}
      
      # permutate current_grp
      # wrong:  permutated_grp <- sample(current_grp)
      
      # compute pvalues etc
      out_resampled <- dplyr::rbind_all(
        lapply(vars, function(v)
          get_row1(y = sample(v), group = current_grp, covariates = covariates, test = test)))
      
      # save in data.frame
      #       results[results$inputvar %in% notnas, colnames(results) %in% columns_results] <- out_resampled
      
      p_rand[rloop] <- min(out_resampled[ , "pvalue"], na.rm = TRUE)
      
      # in last iteration, add corrected pvals
      if(rloop == repetitions)
        R$p_resampling <- sapply(R[ , "pvalue"], function(x) sum(p_rand < x) / repetitions)
      
    }
    
    
  } # here resampling-loop over
  
  # add sampling p-values to dataframe
  R$p_fdr <- p.adjust(R[ , "pvalue"], method = "fdr")
  R$p_bonf <- p.adjust(R[ , "pvalue"], method = "bonferroni")
  
  plotlist <- NULL
  #boxplots with influence of covariates partialed out
  if(no_plots == FALSE && test == "linreg") {
    
    # use only columns that have values and which are below the cutoff
    cols <- which(!is.na(R$pvalue))
    cols <- cols[which((R$pvalue[cols] < boxplot_args$p_cutoff))]
    
    if(length(cols) >= 1) {  # before if(any(!is.na(R$pvalue))) {
      
      # plot and save pdf
      pdf(file = paste(compname, "_", compare[1], "_", compare[2], "_",
                       "boxplots", if(vloop != 0){paste("_", colnames(vargroups)[vloop])},
                       ".pdf", sep=""), onefile=T)
      
      plotlist <- make_boxplots.pairwise(results = R[cols, ], vars = vars[, cols, drop = FALSE], groups = current_grp, 
                             covariates = covariates, boxplot_args = boxplot_args)
      dev.off()
      
      
      # or tell why it failed
    } else if(no_plots == FALSE && verbose) {
      
      notna <- sum(!is.na(R$pvalue))
      below_cutoff <- sum(R$pvalue < boxplot_args$p_cutoff, na.rm = TRUE)
      
      if(notna == 0)
        message(paste0("All p-values NA. Omitting boxplots for ", compare[1], "_vs_", compare[2]))
      else if(below_cutoff == 0)
        message(paste0("No p-values below p_cutoff. Omitting boxplots for ", compare[1], "_vs_", compare[2]))
    }
    
    
    # print volc to pdf
    if((length(R$beta) > 1) && any(!is.na(R$beta))) {
      pdf(paste(compname,"_", compare[1], "_", compare[2], "_", "volcano",
                if(vloop!=0){paste("_", colnames(vargroups)[vloop])}, ".pdf", sep = ""))
      
      plotlist <- append(plotlist, print_volc(datframe = R, volc_args = volc_args))
      dev.off()
      
      # or tell why it failed
    } else if(verbose) {
      
      # below_cutoff <- sum(R$pvalue < p_cutoff, na.rm = TRUE)
      notna <- sum(!is.na(R$pvalue))
      
      if(notna == 0)
        message(paste0("All p-values NA. Omitting volcano plot for ", compare[1], "_", compare[2]))
      else #if(below_cutoff == 0)
        message(paste0("Omitting volcano plot for ", compare[1], "_", compare[2], " for unknown reason. Please contact maintainer."))
    }
  }
  
  # list(R = R, plotlist = plotlist)
  R
}
