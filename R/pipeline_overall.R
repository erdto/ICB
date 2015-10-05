

#' Group comparisons of metabolomics data.
#'
#' Computes Overall F-tests (equivalent) or fisher-tests, for all groups in all metabolites given in "data".
#' 
#' @param data  [data.frame] with metabolites as columns.
#' @param groups  [vector] grouping variable with length matching the number of rows of \code{data}
#' @param compname  [character] prefix for name of output files (excel file, boxplots, volcano plots).
#' @param vargroups [data.frame] with columns being grouping variables of the metabolites in \code{data}. 
#' Used for computing aggregated z-scores for the groups to be compared instead of comparing the single metabolites.
#' @param test [character] what test is to be performed? Linear regression F-test (if no covariates provided, equivalent to Anova) or
#' fishers exact test for comparing the number of NA (below threshold concentration) vs. non-NA values per metabolite.
#' @param repetitions [numeric] Number of repetitions for computing resampling-corrected p-values. Usually, at least 200 are recommended,
#' though higher accuracy can only obtain with greater numbers of repetitions.
#' @param ratios [logical] Should ratios of the metabolites be computed an tested? Will quickly lead to long runtimes.
#' @param ratioop [character] Should ratios be computed by "subtract"ing the metabolite values (for logarithmic data) or by "divide"ing?
#' @return Returns NULL, all output (excel file, pdf with boxplots and volcano plots) is saved into the working directory.
#' @export
pipeline_overall<-  function(data, 
                             groups,
                             refcat, #reference category (control group)
                             compname, 
                             vargroups = data.frame(), 
                             test = c("linreg", "fisher"), 
                             covars = data.frame(), 
                             ratio_args = set_ratio_args(),
                             repetitions = 0,
                             boxplot_args = set_boxplot_args2()) {
  
  covariates <- covars
  # check inputs
  checkmate::assertDataFrame(data, types = "numeric", all.missing = TRUE, 
                             min.rows = 8, min.cols = 1)
  checkmate::assertFactor(groups, empty.levels.ok = FALSE, all.missing = FALSE, len = nrow(data))
  checkmate::assertDataFrame(vargroups, types = c("character", "factor"))
  checkmate::assertCharacter(compname, min.chars = 1, len = 1)
  checkmate::assertDataFrame(covariates, types = c("numeric", "character", "factor"),
                             all.missing = FALSE)
  #   checkmate::assertIntegerish(ncores, lower = 0, len = 1)
  checkmate::assertIntegerish(repetitions, len = 1)
  checkmate::assertChoice(boxplot_args$no_plots, c(TRUE, FALSE))
  
  if(any(table(groups) < 3)) {
    lvltodrop <- names(table(groups)[which(table(groups) <3)])
    
    # leave out rows from this group
    data <- data[groups != lvltodrop , ]
    covariates <- covariates[groups != lvltodrop, ]
    groups <- factor(groups, levels = unique(groups)[!unique(groups) %in% lvltodrop])
    groups <- na.exclude(groups)
    warning(sprintf("Found factor level in 'groups' with less than 3 cases, this will likely not be useful. Found for %s", 
                 names(table(groups)[which(table(groups) <3)])))
  }
  if(missing(refcat))
    refcat <- levels(groups)[1]
  if(!refcat %in% unique(groups))
    stop(sprintf("Argument refcat has to be one of %s", paste0(unique(groups), collapse = ",") ))
  
  ####################
  ####################
  
  vars <- data
  vncols <- NCOL(vargroups)
  ratio_args
  ratios <- ratio_args$ratios
  test <- match.arg(test)

  #create file to save workbook to
  file <- paste(getwd(), "/", compname, ".xls", sep="")
  wb <- XLConnect::loadWorkbook(file, create=TRUE)
  
  ################## vloop 0:vncols #######################
  for(vloop in 0:vncols) {
    
    if(vloop > 0)
      vars <- use_agg_z_df(df = data, vargrp = vargroups[ , vloop])
    
    if(ratios)
      vars <- add_ratios_to_df(df = vars, names = names(vars))#, ratioop = ratio_args$ratio_op)
    
    
    ############# loop over reps #################
    for(rloop in 0:repetitions) {
      if(rloop==1){
        cat("\n", "resampling...", "\n")
      }
      
      
      if(rloop > 0) {
        #permutate y
      }
      
      # loop over metabolites
      results <- data.frame(inputvar = names(vars), dplyr::rbind_all(
        lapply(vars, function(var) {
          get_row3(var, groups, covars, refcat = refcat, test = test)
        })))
      
      
      if(rloop == 0) {
        
        R <- results
        R$p_fdr <- p.adjust(R$p_overall, method="fdr")
        R$p_bonf <- p.adjust(R$p_overall, method="bonferroni")
        
      } else if(rloop != 0) {
        
          p_rand[rloop] <- min(results[ , "p_overall"], na.rm = TRUE)
        
      }
    } # here resampling-loop over   
    
    # add sampling p-values to dataframe
    if(repetitions > 0)
      R$p_sampling <- sapply(results[,"p_overall"], 
                             function(x) sum(p_rand<x)/repetitions)
    
    ################ boxplots ####################
    #boxplots with influence of covars partialed out
    if(boxplot_args$no_plots == FALSE && test=="linreg" && rloop == 0) {
    
      # sort boxplots after pvalues
      if(boxplot_args$sort_boxplots) {
        R <- R[order(R$p_overall, decreasing = FALSE), ]
        vars <- vars[ , as.character(R$inputvar)]
      }
      
      # use only columns that have values and which are below the cutoff
      cols <- which(!is.na(R$p_overall))
      cols <- cols[which((R$p_overall[cols] < boxplot_args$p_cutoff))]
      
      if(length(cols) >= 1) {  # before if(any(!is.na(R$pvalue))) {
        
        # open boxplot pdf
        pdf(file=paste(compname,"_", "multigroup", "_", "boxplots", 
                       if(vloop!=0){paste("_", colnames(vargroups)[vloop])},
                       ".pdf", sep=""), onefile=T)
        
        
        # boxplot_args
        sapply(names(vars)[cols], function(j) {
          if(!is.na(R[R$inputvar == j, "p_overall"][1]))
            print_boxplot2(data.frame(y = vars[, j],
                                      group = groups),
                           name = j,
                           p_overall =  R[R$inputvar == j, "p_overall"][1], 
                           boxplot_args = boxplot_args)
        })
        dev.off()
      }
    }
    # vloop=0 -> metabolite loop
    # add pgains
    
    if(ratios & vloop == 0) {
      R <- compute_pgain(R)
      
      #       R[,"pgain"] <- NA
      #       ratioindices <- (length(varnames)+1):nrow(R)
      #       pratio <- R[ratioindices, "pvalue"]
      #       
      #       p1 <- p2 <- NA
      #       p1 <- sapply(ratio[1, ], FUN = function(x) R[R$inputvar==x, "p_overall"]  )
      #       p2 <- sapply(ratio[2, ], FUN = function(x) R[R$inputvar==x, "p_overall"]  )
      #       
      #       pg <- pratio/pmin(p1,p2)
      #       R[ratioindices, "pgain"] <- pg
    }
    
    # write results to excel file
    SN <- paste(if(vloop==0){"Metabolites"}else{paste(colnames(vargroups)[vloop])}, sep="")
    sheet <- XLConnect::createSheet(wb, name = SN)
    
    # write out, sorted by p_overall
    R <- R[with(R, order(p_overall)), ]
    
    if(test == "fisher") {
      ind <- which(names(R) == "p_overall")
      names(R)[ind] <- "pvalue"
    }
    XLConnect::writeWorksheet(wb, R, sheet = SN)
  } # vloop end
  XLConnect::saveWorkbook(wb, file=file)
  cat("\n", "done.", "\n")
  invisible(0)
}
