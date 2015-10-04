#library(xlsx)
#install.packages("ggplot2")
#install.packages("XLConnect")
# library(ggplot2)
# library(XLConnect)

# load utility functions
# source(file.path(getwd(), "basic_pipelines", "utilities_pipeline.R"))

#' Pipeline for overall F-tests of multiple groups.
#' 
#' 
pipeline_overall2 <- function(data, 
                              groups,
                              refcat, #reference category (control group)
                              compname, 
                              vargroups=data.frame(), 
                              test="linreg", 
                              covars=data.frame(), 
                              ratios=FALSE, 
                              ratioop="subtract", 
                              ratio_cutoff=.1, 
                              repetitions=0,
                              boxplot_args = set_boxplot_args2() ) {
  
  
  vars <- data
  varnames <-  names(data)
  covarnames <- names(covars)
  ngrp <- length(unique(groups)) # number of grps
  group <- unique(groups)
  group[1] <- refcat           #reordering groups
  group[2:ngrp] <- group[group!=refcat]
  vargroups <- na.exclude(vargroups)
  vncols <- ncol(vargroups)
  
  if(ratios)
    data <- add_ratios_to_df(df = data, names = varnames)
  
  for(vloop in 0:vncols) {
    
    if(vloop != 0){
      cat("\r", "vargroup-loop:", vloop, "out of", vncols) # , "\n")
      
      # y != 0 means we're in a vargroup-loop, new data matrix with agg z scores will be used
      vars <- use_agg_z_df(df = data, vargrp = vargroups[ , vloop])
      
      if(ratios){
        varnames <- names(vars)
        vars <- add_ratios_to_df(df = vars, names = varnames)
      }
    }
    
  
    #create file to save workbook to
    file <- paste(getwd(), "/", compname, ".xls", sep="")
    # alt: 
    #wb <- createWorkbook(type="xls")
    wb <- XLConnect::loadWorkbook(file, create=TRUE)
    
    
    
    # open boxplot pdf
    pdf(file=paste(compname,"_", "multigroup", "_", "boxplots", if(vloop!=0){paste("_", colnames(vargroups)[vloop])},".pdf", sep=""), onefile=T)
    
    # resampling: vector for correction
    p_rand <- matrix(NA)
    
    for(rloop in 0:repetitions) {
      if(rloop==1){
        cat("\n", "resampling...", "\n")
      }
      #         if((rloop/(repetitions)*100)%%10==0){
      #           cat( "  ", round((rloop/(repetitions)*100), digits=1), "%", "  ", sep="")
        #         }
        
        #initialize matrix results where results will be saved and which is written into excel sheets later
        results <- as.data.frame(matrix(NA, nrow=ncol(vars),ncol= 6+3*ngrp))
        names(results) <- c("inputvar", "N", paste("n", 1:ngrp, sep=""), 
                            "p_Overall", paste("beta", group, sep="_"), 
                            paste("pvalue", group, sep="_"), "p_fdr", "p_bonf", "p_sampling") 
        results$inputvar <- colnames(vars)
        #results_list <- (dim(c(dim(results),ncol(grpcomb))))
        
        
        
        for(j in 1:length(vars)) {
          if(rloop==0) {
#             if((j/length(vars)%%10)==0){cat("  ",round(((j)/(length(vars)))*100, digits=1), "%", "  ", sep="")
#             } 
            
            if(nrow(covars) > 0)
              X <- data.frame(vars[,j], groups, covars)
            else
              X <- data.frame(vars[,j], groups)
            X <- within(X, groups <- relevel(groups, ref = refcat))
            names(X) <- c("y","group", covarnames)
            
            if(test=="fisher"){
              Y <- matrix(c(sapply(group, function(x) sum(is.na(X[X$group==x,]))), sapply(group, function(x) sum(!is.na(X[X$group==x,]))))
                          , byrow=T, ncol=length(group), nrow=2, dimnames = list(missing = c("NA", "non-NA")))
            }
          
          #X2 <- data.frame(vars[,j], rep(NA, times=length(vars[,j])), covars) # oder as.matrix(covars) ?
          #str(X2)
          X <- X[which(X$group!="NA"),]
          
          #abfragen, ob f?r jede gruppe mindestens 2 cases nicht NA sind:
          sum_cases <- NULL
          sum_cases <- sapply(1:ngrp, function(x) { nrow(X[X$group == group[x] & !is.na(X$y),]) })
#           if(any)
#           sum_cases_covars <- NULL
#           sum_cases_covars <- sapply(1:ngrp, function(x) { nrow(X[X$group == group[x] & !is.na(X$y),]) })
#         
          


          if(test=="linreg"){
            err <- tryCatch(
              lm(y ~ . , data=X)
          , error = function(e){ TRUE })}
          
          else if(test=="fisher"){
            err <- tryCatch( 
              ft <- fisher.test(Y)
          , error = function(e){ TRUE })
          }

          if(class(err)!= "logical"){err <- FALSE}

          
          # if there are less that 2 data points, NA:
          if(any(sum_cases < 3) || err==TRUE){
            next
          }
          #if resampling, mix up groups, then lm :
          else if(rloop!=0) { # && ((nrow(X[X$G==0 & !is.na(X$y),]) > 1) & (nrow(X[X$G==1 & !is.na(X$y),]) > 1))){
            #mix it up:
            for(p in 1:ngrp){
              X[X$group==group[p],"group"] <- sample.int(n=ngrp, size=length(X[X$group==unique(groups)[p],"group"]), replace=T)-1
            }
            test.out <- lm(y ~ . , data=X)
            results[results$inputvar==colnames(vars)[j], "pvalue"] <- summary(test.out)$coef[2,4]
          }
          else{ 
            # here,we are before the resampling and 
            #the data is prepped for regression and plotting
            
            if(test=="linreg"){
              test.out <- lm(y ~ . , data=X)
              p <- summary(test.out)$coef[1:ngrp,4]
              
              f <- summary(test.out)$fstatistic
              p_f <- pf(f[1],f[2],f[3],lower.tail=F)
              
              beta <- test.out$coef[1:ngrp]
              CIlower <- confint(test.out)[1:ngrp,1]
              CIupper <- confint(test.out)[1:ngrp,2]
              
              # print results to df:
              results[results$inputvar==colnames(vars)[j],]$N <- nrow(X)
              results[results$inputvar==colnames(vars)[j],]$n1 <- sum(X$group==refcat)
              results[results$inputvar==colnames(vars)[j], paste("beta", refcat ,sep="_")] <- unname(beta[1])
              results[results$inputvar==colnames(vars)[j], "p_Overall"] <- p_f
              results[results$inputvar==colnames(vars)[j], paste("pvalue", "_", refcat, sep="")] <- p[1]
              
              
              for(n in 2:ngrp){
                results[results$inputvar==colnames(vars)[j], paste("beta", group[n] ,sep="_")] <- unname(beta[n])
                results[results$inputvar==colnames(vars)[j],paste("CI_lower", refcat, sep="_")] <- unname(CIlower[1])
                results[results$inputvar==colnames(vars)[j],paste("CI_upper", refcat ,sep="_")] <- unname(CIupper[1])
                results[results$inputvar==colnames(vars)[j],paste("CI_lower", group[n] ,sep="_")] <- unname(CIlower[n])
                results[results$inputvar==colnames(vars)[j],paste("CI_upper", group[n] ,sep="_")] <- unname(CIupper[n])
                
                results[results$inputvar==colnames(vars)[j], paste("n",n, sep="")] <-  sum(X$group==levels(X$group)[n])
                results[results$inputvar==colnames(vars)[j], paste("pvalue", group[n], sep="_")] <- p[n]
              }
            }
            else if(test=="fisher"){
              ft <- fisher.test(Y, simulate.p.value = TRUE, B = 1e5)
              
              # print results to df:
              results[results$inputvar==colnames(vars)[j],]$N <- nrow(X)
              results[results$inputvar==colnames(vars)[j], beta] <- NA
              results[results$inputvar==colnames(vars)[j], pvalue] <- ft$p.value
              results[results$inputvar==colnames(vars)[j],]$CIlower <- ft$conf.int[1]
              results[results$inputvar==colnames(vars)[j],]$CIupper <- ft$conf.int[2]
            }
            
            #boxplots with influence of covars partialed out
            
            
            if(test=="linreg" && rloop==0 && (ratios==FALSE | (ratios==TRUE & p < ratio_cutoff))) {
              
              boxplot_args
              boxplot_theme = boxplot_args$boxplot_theme
              boxplot_colors = boxplot_args$boxplot_colors
              boxplot_fill = boxplot_args$boxplot_fill
              boxplot_size = boxplot_args$boxplot_size
              boxplot_shape = boxplot_args$boxplot_shape
              boxplot_alpha = boxplot_args$boxplot_alpha

              
              title <- sprintf("Overall test p = %.6g",p_f)
              X <- na.exclude(X) 
              
              #for(x in 1:ncol(covars)){if(length(levels(X[,x]))>5){X[,x] <<-as.numeric(as.character(X[,x]))}}
              X$yb <- lm(y ~ . -group, data=X)$res
              #summary(lm(y ~ . -G, data=X))
              
              theme <- theme_set(theme_minimal())
              b <- qplot(factor(group), yb, data=X, geom="boxplot", colour=factor(group)) + 
                ylab(colnames(vars)[j]) +  ggtitle(title) +
                geom_jitter(aes(color=factor(group))) + 
                scale_x_discrete("Group", labels= group) + 
                theme(legend.position="none") + 
                scale_colour_manual(values=boxplot_colors) +
                scale_size_manual(values= boxplot_size)  +
                scale_alpha_manual(values = boxplot_alpha) +
                scale_shape_manual(values = boxplot_shape) +
                scale_fill_manual(values = boxplot_fill)
              
              if(class(boxplot_theme)[1]=="theme"){b <- b + boxplot_theme}  #+ boxplot_layers #theme_update(theme, unlist(boxplot_theme))
              
              print(b)
            }
          }
        }} # vars loop over
          
        if(rloop == 0){
          R <- results
          dev.off()
        }
        else if(rloop != 0){
          p_rand[rloop] <- min(results[,"pvalue"], na.rm=T)
        }
        
      } # here resampling-loop over    
      
      # add sampling p-values to dataframe
      R$p_fdr <- sapply(R$p_Overall, function(x) p.adjust(x, method="fdr"))
      R$p_bonf <-sapply(R$p_Overall, function(x) p.adjust(x, method="bonferroni"))
      R$p_sampling <- sapply(results[,"p_overall"], function(x) sum(p_rand<x)/repetitions)
      
      
      
      # vloop=0 -> metabolite loop
      
      if(ratios==TRUE & vloop == 0) {
        R <- compute_pgain(R)
#         R[,"pgain"] <- NA
#         ratioindices <- (length(varnames)+1):nrow(R)
#         pratio <- R[ratioindices, "pvalue"]
#         
#         p1 <- NA
#         p1 <- sapply(ratio[1, ], FUN = function(x) R[R$inputvar==x, "p_Overall"]  )
#         p2 <- NA
#         p2 <- sapply(ratio[2, ], FUN = function(x) R[R$inputvar==x, "p_Overall"]  )
#         
#         pg <- pratio/pmin(p1,p2)
#         R[ratioindices,"pgain"] <- pg
      }
      
      # not rly needed
      #results_list[[i]] <- R 
      
      # write results to excel file
      SN <- paste(if(vloop==0){"Metabolites"}else{paste(colnames(vargroups)[vloop])}, sep="")
      sheet <- XLConnect::createSheet(wb, name = SN)
      
      #alt:
      #addDataFrame(results, sheet=sheet, col.names=TRUE, row.names=F)
      XLConnect::writeWorksheet(wb, R, sheet = SN)
  }
  XLConnect::saveWorkbook(wb, file=file)
  cat("\n", "done.")
}

