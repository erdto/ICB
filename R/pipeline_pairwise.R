# options(java.parameters = "-Xmx1024m")
# library(ggplot2)
# library(XLConnect)
# library(parallel)
# library(checkmate)

#' Pairwise group comparisons for metabolomics data.
#'
#' Computes pairwise t-tests (or linear regressions when covariates present) or fisher-tests, of the difference of all groups in all metabolites given in "data".
#' This function is intended for exploratory analyses. The type-I-error-rate can be controlled with resampling or bonferroni-/FDR-corrections. 
#' 
#' @param data  [data.frame] with metabolites as columns.
#' @param groups  [vector] grouping variable with length matching the number of rows of \code{data}
#' @param covars [data.frame] with variables to be controlled for.
#' @param compname  [character] prefix for name of output files (excel file, boxplots, volcano plots).
#' @param vargroups [data.frame] with columns being grouping variables of the metabolites in \code{data}. 
#' Used for computing aggregated z-scores for the groups to be compared instead of comparing the single metabolites.
#' @param test [character] what test is to be performed? Linear regression (if no covariates provided, equivalent to t-test) or
#' fishers exact test for comparing the number of NA (below threshold concentration) vs. non-NA values per metabolite.
#' @param repetitions [integer] Number of repetitions for computing resampling-corrected p-values. Usually, at least 200 are recommended,
#' though higher accuracy can only obtain with greater numbers of repetitions.
#' @param ratio_args see \link{\code{set_ratio_args}}
#' @param boxplot_args see \link{\code{set_boxplot_args}}
#' @param volc_args see \link{\code{set_volc_args}}
#' @param no_plots [logical] Should all plots be omitted? This will speed up the function considerably.
#' @param verbose [logical] Should progress be printed?
#' @param ncores [integer] number of CPUs to be used. If >1, computations will be parallelized.
#' @return Returns 0, all output is saved in an excel file and plots in pdfs in the working directory.
#' @export
pipeline_pairwise <- function(data, groups, compname,
                              vargroups = data.frame(),
                              covars = data.frame(),
                              test = c("linreg", "fisher"),
                              ratio_args = set_ratio_args(),
                              repetitions = 0,
                              boxplot_args = set_boxplot_args(),
                              volc_args = set_volc_args(),
                              no_plots = FALSE,
                              verbose = TRUE,
                              ncores = 1) {
  compare = NULL
  ratios <- ratio_args$ratios
  ratioop <- ratio_args$ratio_op 
  
  checkmate::assertDataFrame(data, types = "numeric", all.missing = TRUE)
  if(nrow(covars) != 0)
    checkmate::assertDataFrame(covars, types = c("numeric", "factor"), all.missing = FALSE)
  checkmate::assertFactor(groups, len = nrow(data), all.missing = FALSE)
  
  if((!class(groups) %in% c("character", "factor")) || length(groups) != NROW(data))
    stop("groups has to of type factor and of length nrow(data).")
  
  if(NROW(data) == 0 || NROW(groups) == 0) {
    warning("Empty input in 'data' or 'groups'. Returning NA.")
    return(NA)
  }
  if(!(NROW(data) == NROW(groups)))
    stop("Unequal number of samples in 'data' and 'groups'.")
  if(!(NROW(covars) == 0 | NROW(data) == NROW(covars)))
    stop("Unequal number of rows in 'data' and 'covars'.")
  if(!(NROW(covars) == 0 | NROW(covars) == NROW(groups)))
    stop("Unequal number of rows in 'groups' and 'covars'.")
  
  
  data <- as.data.frame(data)
  covariates <- as.data.frame(covars)
  vargroups <- as.data.frame(vargroups)
  stopifnot(is.data.frame(data) && is.data.frame(vargroups)
            && is.data.frame(covariates))
  if(ratios == TRUE & no_plots == FALSE)
    warning(paste("Warning, ratios = TRUE and no_plots = FALSE, this could take a while!", "\n",
                  "Using p_cutoff argument is recommended when ratios is on."))
  
  varnames <-  colnames(data)
  covarnames <- colnames(covariates)
  #   vargroups <- na.exclude(vargroups)
  vncols <- NCOL(vargroups)
  test <- match.arg(test)  
  # if not specified in compare-argument, do comparisons on all
  # pairwise combinations of the groups
  #   if(is.null(compare)) {
  compare <- combn(levels(groups), m = 2, simplify = FALSE)
    
  
  #   }
  
  # We use 'data' inside the loops later on, so we save the original in
  # 'data1' to call 'use_agg_z_df' with at the beginning of the vargroup-loops
  data1 <- data
  
  if(ratios==TRUE){
    data <- add_ratios_to_df(df = data, names = varnames, ratioop = ratioop)
  }
  
  # Prepare output excel-file
  file <- paste(getwd(), "/", compname, "_pipeline",
                if(test == "fisher"){"_fisher"}, ".xls", sep="")
  wb <- XLConnect::loadWorkbook(file, create=TRUE) # xls-alternative: createWorkbook(type="xls")
  sheetnamez <- NULL
  
  # Specify to NOT wrap the text
  cs <- XLConnect::createCellStyle(wb)
  XLConnect::setWrapText(cs, wrap = FALSE)
  
  for(vloop in 0:vncols) {
    
    if(vloop != 0) {
      cat("\r", "vargroup-loop:", vloop, "out of", vncols) # , "\n")
      if(vloop == vncols) cat("\n")
      
      # y != 0 means we're in a vargroup-loop, new data matrix with agg z scores will be used
      data <- use_agg_z_df(df = data1, vargrp = vargroups[ , vloop])
      
      if(ratios==TRUE)
        data <- add_ratios_to_df(df = data, names = names(data))
    }
    
    # Compare groups does everything:
    
    if(.Platform$OS.type != "windows") {
      results_list <- parallel::mclapply(compare, function(x) 
        compare_grps(df = data, groups = factor(groups, levels = x), covariates = covariates, compname = compname,
                     vloop = vloop, vargroups = vargroups,
                     repetitions = repetitions, test = test,
                     boxplot_args = boxplot_args, volc_args = volc_args, no_plots = no_plots,
                     verbose = verbose),
        mc.cores = ncores)
      
    } else if(ncores == 1) {
      
      results_list <- lapply(compare, function(x) 
        compare_grps(df = data, groups = factor(groups, levels = x), covariates = covariates, compname = compname,
                     vloop = vloop, vargroups = vargroups,
                     repetitions = repetitions, test = test,
                     boxplot_args = boxplot_args, volc_args = volc_args, no_plots = no_plots,
                     verbose = verbose))
      
      
    } else {
      
      cl <- parallel::makePSOCKcluster(rep("localhost", ncores))
      on.exit(parallel::stopCluster(cl))
      
      # functions called by compare_grps have to be exported explicitly
      parallel::clusterExport(cl = cl, varlist = c("compare_grps", "set_boxplot_args", "set_volc_args", 
                                                   "make_boxplots.pairwise", "print_boxplot", "print_volc", "get_row1"),  envir = environment())
      
      # if your function does not call "require(ggplot2)"
      parallel::clusterEvalQ(cl, library(ggplot2))
      
      results_list <- parallel::parLapply(cl = cl, X = compare, function(X) 
        compare_grps(df = data, groups = factor(groups, levels = X), covariates = covariates,
                     vloop = vloop, vargroups = vargroups, compname = compname,
                     repetitions = repetitions, test = test,
                     boxplot_args = boxplot_args, volc_args = volc_args, no_plots = no_plots,
                     verbose = verbose))
    }
   
#     plotlist <-  lapply(results_list, function(i) i$plotlist), nm = compare)
#     results_list <- lapply(results_list, function(i) i$R)

    for(page in seq_along(compare)) {
      
      # if ratios are used
      # split varnames into two columns and compute p_gain
      if(ratios) {
        ssplit <- strsplit(as.character(results_list[[page]][, "inputvar"]), "/")
        ratioind <- (length(varnames)+1):nrow(results_list[[page]])
        
        results_list[[page]][, "inputvar"] <- sapply(ssplit, function(x) x[[1]]) 
        results_list[[page]][ratioind, "denominator"] <- 
          unlist(lapply(ssplit[ratioind], function(x) ifelse(length(x)>1,x[[2]], NA)))
        
        # reorder cols
        results_list[[page]] <- results_list[[page]][, 
                                                     c("inputvar", "denominator", 
                                                       names(results_list[[page]])[!names(results_list[[page]]) %in% 
                                                                                     c("inputvar", "denominator")])
                                                     ]  
        # add "pgain" column to output
        results_list[[page]][, "pgain"] <- compute_pgain(results_list[[page]])
        
      }
      
      # print out xls sheets
      SN <- paste(compare[[page]][1], "_", "vs.","_", compare[[page]][2],
                  if(vloop != 0){paste("_", colnames(vargroups)[vloop])}, sep = "")
      
      # save the sheet names of the original comparisons, need them for reordering
      sheetnamez <- append(sheetnamez, SN)
      
      if(nchar(SN) > 30) {
        stopifnot(SN == sheetnamez[page + (length(compare) * vloop)])
        warning(
          sprintf("Name of sheet %d not allowed to be longer than 30 characters. Cutting off some.", 
                  page + (length(compare) * vloop))
          
        )
        if(vloop > 0) {
          # only keep vargroup label, as sheets are later ordered after the comparisons
          # i.e. A_vs._B, then A_vs.B_subgroup before, A_vs._C etc.
          SN <- paste0(colnames(vargroups)[vloop], "_", page)
          sheetnamez[page + (length(compare) * vloop)] <- SN
        } else {
          stop(sprintf("Group names in comparison: %s too long for print out in excel sheet.", SN))
        }
      }
      sheet <- XLConnect::createSheet(wb, name = SN)
      
      # save data sorted by pvalue
      XLConnect::writeWorksheet(wb, results_list[[page]][with(results_list[[page]], order(pvalue)), ], 
                                sheet = SN)
      # set cell style
      XLConnect::setCellStyle(wb, sheet = SN, row = 2:(2+length(results_list[[page]][["inputvar"]])), col = 1,
                              cellstyle = cs)
      
    }
  } # vloop over
  
  # if there are "pathway"-sheets, order them
  if(vncols > 0) {
    stopifnot((length(sheetnamez) == length(compare) * (vncols+1)))
    # sheetnamez_0
    #     paste(compare[[page]][1], "_", "vs.","_", compare[[page]][2],
    #           if(vloop != 0){paste("_", colnames(vargroups)[vloop])}, sep = "")
    #     sheetz <- sheetnamez 
    # XLConnect::getSheets(wb)
    # count for sheet sequence
    count <- 1
   # sequence of sheet names
    index_mat <- 
    matrix(1:length(compare), byrow = TRUE, ncol = length(compare), nrow = vncols+1) + 
      matrix(c(0, length(compare) * 1:vncols), nrow = vncols+1, ncol = length(compare))
    
    # before sapply(1:length(compare), function(i)
    # c(i, i + length(compare) * 1:vncols))
       
    # set sheet positions
    sapply(as.numeric(index_mat), function(j) {
      XLConnect::setSheetPos(object = wb, sheet = rev(sheetnamez)[j], pos = count)
      count <- count + 1
      NULL
    })
  }
  
  # save workbook to disk
  XLConnect::saveWorkbook(wb, file=file)
  cat("\n", "done.", "\n")
  
#   list(results_list = results_list,
#        plotlist = plotlist)
  # results_list
  invisible(0)
}



