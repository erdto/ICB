
#' Utility function for performing the tests
#' 
#' Inputs:  y response [numeric vector] , 
#'          group grouping [factor vector] ,
#'          covariates [data.frame] with columns being covariates
#'          test [choice] which test to perform
#'          
#' Outputs:  Row with sample sizes, coefficient, p-value, CI. 
get_row1 <- function(y, groups, covariates = data.frame(), test = c("linreg", "fisher")) {
  stopifnot(is.data.frame(covariates))
  stopifnot(!is.null(y))
  group <- groups
  
  test <- match.arg(test)
  
  number_samples <- sum(!is.na(y))
  grp1 <- levels(group)[!is.na(levels(group))][1]
  grp2 <- levels(group)[!is.na(levels(group))][2]
  
  
  if(nrow(covariates) > 0) {
    # remove Nas from covars
    Nas <- union(which(is.na(y)), which(is.na(group)))
    covariates2check <- covariates[!1:NROW(covariates) %in% Nas, ]
    
    # check only after NAs have been removed!
    for(l in names(covariates2check)) {
#       if(length(covariates2check[, l]) == length(unique(covariates2check[ , l])) && is.factor(covariates2check[,l]))
#         stop(sprintf("Covariate %s is encoded as factor with all unique values. Cannot estimate anything from this.", l))
      if(nrow(covariates2check) > 4 & length(unique(covariates2check[, l])) == 1) {
        warning(sprintf("Covariable %s has only one unique value. Skipping this test." , l))
        
        df <- as.data.frame(y = y, group = group, covariates)
#         df$y <- y
        n1 <- nrow(na.exclude(df[group == grp1, ]))
        n2 <- nrow(na.exclude(df[group == grp2, ]))
        
        return(data.frame(
          N = number_samples, 
          n1 = n1, 
          n2 = n2, 
          beta = NA, 
          pvalue = NA, 
          CI_lower = NA, 
          CI_upper = NA))
      }
    }
  }
  
  if(nrow(covariates)!=0 && nrow(covariates) != length(y)) {
    print(y)
    print(covariates)
    stop("number of rows in covariates and response differs.")
  }
  
  
  if(test == "fisher") {
    #     number_samples <- length(y)
    #     missing <- sum(is.na(y))
    
    df <- data.frame(na = is.na(y), group = group)
    # remove rows from NA groups
    df <- df[which(!is.na(df[,"group"])), ]
    
    n1 <- sum(df$group == grp1 & !df$na)
    n2 <- sum(df$group == grp2 & !df$na)
    
    na_1 <- sum(df[which(df$group == grp1), "na"])
    na_2 <- sum(df[which(df$group == grp2), "na"])
    N <- n1 + n2 + na_1 + na_2
    #     N_na <- n1_na + n2_na
    expected_n1 <- ((n1 + na_1) / N) * ((n1 + n2) / N) * N
    expected_n2 <- ((n2 + na_2) / N) * ((n1 + n2) / N) * N 
    
    
    #     n1 <- sum(!is.na(y[group == grp1]))
    #     n2 <- sum(!is.na(y[group == grp2]))
    #     n1_na <- sum(is.na(y[which(group == grp1)]))
    #     n2_na <- sum(is.na(y[which(group == grp2)]))
    
    # remove NAs from groups
    fishy <- try(fisher.test(as.factor(is.na(y[!is.na(group)])), as.factor(group[!is.na(group)]), conf.int = TRUE), silent = TRUE)
    
    # if error return empty line
    if(inherits(fishy, "try-error"))
      return(
        data.frame(n1 = n1, 
                   na_1 = na_1, 
                   n2 = n2, 
                   na_2 = na_2,
                   exp_n1 = expected_n1,
                   exp_n2 = expected_n2,
                   pvalue = NA, 
                   oddsratio = NA,
                   CI_lower = NA,
                   CI_upper = NA, row.names = NULL))
    
    if(is.null(fishy$conf.int) || any(is.na(fishy$conf.int)))
      CI_lower <- CI_upper <- NA
    else {
      CI_lower <- fishy$conf.int[1]
      CI_upper <- fishy$conf.int[2]
    }
    if(!is.null(fishy$estimate))
      beta <- fishy$estimate
    else
      beta <- NA
    
    return(data.frame(
      n1 = n1,
      na_1 = na_1,
      n2 = n2,
      na_2 = na_2, #
      exp_n1 = expected_n1,
      exp_n2 = expected_n2,
      oddsratio = beta,
      pvalue = fishy$p.value,
      CI_lower = CI_lower,
      CI_upper = CI_upper, row.names = NULL))
    
  } else {
    
#     number_samples <- sum(!is.na(y))
#     grp1 <- unique(group)[!is.na(unique(group))][1]
#     grp2 <- unique(group)[!is.na(unique(group))][2]
    
    if(NROW(covariates) == 0) {
      df <- data.frame(group = group, y = y)
      df <- na.exclude(df)
      n1 <- nrow(df[df$group == grp1, ])
      n2 <- nrow(df[df$group == grp2, ])
      
      lmy <- try(lm(y ~ group), silent = TRUE)
      
    } else {
      
      df <- data.frame(y = y, group = group, covariates)
      df <- na.exclude(df)
      n1 <- nrow(df[df$group == grp1, ])
      n2 <- nrow(df[df$group == grp2, ])
      
      lmy <- try(lm(y ~., data = df), silent = TRUE)
    }
  }
  
  if(inherits(lmy, "try-error") || (n1 < 2 || n2 < 2)) {
    return(data.frame(
      N = number_samples, 
      n1 = n1, 
      n2 = n2, 
      beta = NA, 
      pvalue = NA, 
      CI_lower = NA, 
      CI_upper = NA))
  } else {
    sumy <- summary(lmy)
    confy <- confint(lmy)
    return(data.frame(N = number_samples,
                      n1 = n1,
                      n2 = n2,
                      beta = sumy$coef[2, 1],
                      pvalue = sumy$coef[2, 4],
                      CI_lower = confy[2, 1],
                      CI_upper = confy[2, 2], row.names = NULL))
    
  }
}

###############################################################################
#########
#' For pipeline_continuous
get_row2 <- function(X, name, covars, test) {
  stopifnot(is.data.frame(X) && all(c("y", "x") %in% names(X)))
  
  if(nrow(covars) > 0)
    X <- data.frame(y = X$y, name = X$x, covars)
  else
    X <- data.frame(y = X$y, name = X$x)
  names(X) <- c("y", name, names(covars))
  
  if(test=="linreg") {
    if(sum(!is.na(X$y)) < 3)
      return(data.frame(
        predictor = name,
        N = sum(!is.na(X$y)),
        beta = NA,
        pvalue = NA,
        CIlower = NA,
        CIupper = NA, row.names = NULL))
    
    test.out <- lm(y ~ . , data = X)
    
    # just return everything
    sumy <- summary(test.out)
    return(data.frame(
      predictor = name,
      N = sumy$df[1] + sumy$df[2], 
      beta = test.out$coef[name],
      pvalue = sumy$coef[name, 4],
      CIlower = confint(test.out)[name, 1],
      CIupper = confint(test.out)[name, 2]))
    
    
    #     p <- summary(test.out)$coef[prednames[i], 4]
    #     beta <- test.out$coef[prednames[i]]
    #if(is.infinite(-log10(p))){beta <- NA}
    #               results[results$inputvar==colnames(vars)[j],]$CIlower <<- confint(test.out)[prednames[i],1]
    #               results[results$inputvar==colnames(vars)[j],]$CIupper <<- confint(test.out)[prednames[i],2]
  } else if(test == "fisher") {
    X$y <- is.na(X$y)
    
    if(sum(X$y) == 0)
      return(data.frame(
        predictor = name,
        N = nrow(X),
        n_NA = 0,
        beta = NA,
        pvalue = NA,
        CIlower = NA,
        CIupper = NA))
    
    test.out <- glm(y ~ . , data = X, family = "binomial")
    sumy <- summary(test.out)
    beta <- test.out$coef[name]
    return(data.frame(predictor = name,
                      N = sumy$df[1] + sumy$df[2], 
                      n_NA = sum(X$y),
                      beta = exp(beta) / (1 + exp(beta)),
                      pvalue = sumy$coef[name, 4],
                      CIlower = suppressMessages(confint(test.out)[name, 1]),
                      CIupper = suppressMessages(confint(test.out)[name, 2])))
  }
}

###########################################################################################
##############

#' For pipeline_overall
#'
#'  
#' For pipeline_overall
#'
#'  
get_row3 <- function(y, groups, covariates, refcat, test) {
  stopifnot(is.factor(groups))
  
  groups <- relevel(groups, ref = refcat)
  ngrps <- length(levels(groups))
  number_samples <- sum(!is.na(y))
  
  # check number of cases for each grp
  sum_cases <- sapply(levels(groups), function(grp) 
    sum(groups == grp & !is.na(y)))
  
  
  # do test
  if(test == "linreg") {
    
    # input check first:
    # if all NA return only sample size 
    if(any(sum_cases < 3)) {
      # # print results to list:
      rval <- data.frame(N = sum(sum_cases),
                         p_overall = NA,
                         n1 = sum_cases[1])
      group = levels(groups)
      names(rval)[3] <- paste0("n_", group[1])
      
      for(i in 2:ngrps) {
        rval[, paste0("n_", group[i])] <- sum_cases[i]
      }
      return(rval)
    }
    
    
    if(nrow(covariates) > 0)
      X <- data.frame(y = y, group = groups, covariates)
    else
      X <- data.frame(y = y, group = groups)
    
    # relevel
    X <- within(X, group <- relevel(group, ref = refcat))
    
    
    test.out <- summary(lm_out <- lm(y ~ . , data=X, 
                                     contrasts = list(group = contr.treatment)))
    f <- test.out$fstatistic
    p_f <- pf(f[1],f[2],f[3],lower.tail=F)
    
    #     group = lm_out$xlevels$group,
    #     N = sum_cases,
    #     beta = dummy.coef(lm_out)$group,
    #     CI = confint.lm(lm_out)[1:ngrps, ],
    #     pvalue = coef(test.out)[1:ngrps,4],
    #     p_overall = rep(p_f, ngrps)))
    
    
    # # print results to list:
    rval <- data.frame(N = sum(sum_cases),
                       p_overall = p_f,
                       n1 = sum_cases[1])
    group = lm_out$xlevels$group
    names(rval)[3] <- paste0("n_", group[1])
    
    for(i in 2:ngrps) {
      rval[, paste0("n_", group[i])] <- sum_cases[i]
    }
    
    return(rval)
    #       return(data.frame(
    #       group = lm_out$xlevels$group,
    #       N = sum_cases,
    # #       beta = dummy.coef(lm_out)$group,
    # #       CI = confint.lm(lm_out)[1:ngrps, ],
    # #       pvalue = coef(test.out)[1:ngrps,4],
    #       p_overall = rep(p_f, ngrps)))
    
    
  } else if (test == "fisher") {
    
    # X <- data.frame(y = as.factor(as.numeric(is.na(y))), group = groups)
#     
#     if(all(X$y == 0) || all(X$y == 1))
#       return(data.frame(
#         predictor = name,
#         N = nrow(X),
#         n_NA = 0,
#         beta = NA,
#         pvalue = NA,
#         CIlower = NA,
#         CIupper = NA))
# #     
#     test.out <- glm(y ~ ., data = X, family = "binomial")
#     sumy <- summary(test.out)
    
#     beta <- test.out$coef[name]
#     return(data.frame(predictor = name,
#                       N = sumy$df[1] + sumy$df[2], 
#                       n_NA = sum(X$y),
#                       beta = exp(beta) / (1 + exp(beta)),
#                       pvalue = sumy$coef[name, 4],
#                       CIlower = suppressMessages(confint(test.out)[name, 1]),
#                       CIupper = suppressMessages(confint(test.out)[name, 2])))
#     
    X <- data.frame(y = y, group = groups)
    
    Y <- matrix(c(sapply(levels(groups), function(x) sum(is.na(X[X$group == x, "y"]))), 
                  sapply(levels(groups), function(x) sum(!is.na(X[X$group == x, "y"]))))
                , byrow=TRUE, ncol=ngrps, nrow=2, dimnames = list(missing = c("NA", "non-NA")))
    
    ft <- try(fisher.test(Y, simulate.p.value = TRUE, B = 1e5), silent = TRUE)
    p <- if(inherits(ft, "try-error")) {
      NA
    } else {
      ft$p.value
    }
    
    # check number of cases for each grp
    sum_cases <- sapply(levels(groups), function(grp) 
      sum(groups == grp & !is.na(y)))
    sum_nas <- sapply(levels(groups), function(grp) 
      sum(groups == grp & is.na(y)))
    N <- sum(sum_cases) + sum(sum_nas)
    expected_n <- sapply(levels(groups), function(grp)
      ((sum_cases[grp] + sum_nas[grp]) / N) * ((N) / N) * N)
    # TODO
    #     expected_cases <- sapply()
    
    # # print results to list:
    rval <- data.frame(N = N,
                       p_overall = p,
                       n1 = sum_cases[1],
                       n1_exp = expected_n[1],
                       na1 = sum_nas[1])
    group = levels(groups)
    names(rval)[3] <- paste0("n_", group[1])
    names(rval)[4] <- paste0("n_exp_", group[1])
    names(rval)[5] <- paste0("na_", group[1])
    
    for(i in 2:ngrps) {
      rval[, paste0("n_", group[i])] <- sum_cases[i]
      rval[, paste0("n_exp_", group[i])] <- expected_n[i]
      rval[, paste0("na_", group[i])] <- sum_nas[i]
    }
    
    return(rval)
    #     return(data.frame(
    #       variable = unique(groups),
    #       N = sum_cases,
    #       NAs = sum_nas,
    #       pvalue = rep(p, ngrps)))
    
    
  }
}



# get_row3 <- function(y, groups, covariates, refcat, test) {
#   stopifnot(is.factor(groups))
#   
#   ngrps <- length(unique(groups))
#   number_samples <- sum(!is.na(y))
#   
#   # input check first:
#   # if all NA return only sample size 
#   
#   # do test
#   if(test == "linreg") {
#     # check number of cases for each grp
#     sum_cases <- sapply(unique(groups), function(grp) 
#       sum(groups == grp & !is.na(y)))
#     
#     if(any(sum_cases < 3))
#       return(data.frame(
#         group = unique(groups),
#         N = sum_cases
#       ))
#     
#     
#     if(nrow(covariates) > 0)
#       X <- data.frame(y = y, group = groups, covariates)
#     else
#       X <- data.frame(y = y, group = groups)
#     
#     # relevel
#     X <- within(X, group <- relevel(group, ref = refcat))
#     
#     
#     test.out <- summary(lm_out <- lm(y ~ . , data=X, 
#                                      contrasts = list(group = contr.treatment)))
#     f <- test.out$fstatistic
#     p_f <- pf(f[1],f[2],f[3],lower.tail=F)
#     
# #     group = lm_out$xlevels$group,
# #     N = sum_cases,
# #     beta = dummy.coef(lm_out)$group,
# #     CI = confint.lm(lm_out)[1:ngrps, ],
# #     pvalue = coef(test.out)[1:ngrps,4],
# #     p_overall = rep(p_f, ngrps)))
#     sum_cases    
# for(i in 1:ngrps)
#   
#   rownamess <- c("N", "n1", "p_overall", paste("beta", refcat ,sep="_"),)
# # # print results to list:
#     rval <- data.frame(N = sum(sum_cases),
#                  n1 = sum_cases[1],
#                  paste("beta", refcat ,sep="_") = dummy.coef(lm_out)$group[1])
# #     results[results$inputvar==colnames(vars)[j],]$N <- nrow(X)
# #     results[results$inputvar==colnames(vars)[j],]$n1 <- sum(X$group==refcat)
# #     results[results$inputvar==colnames(vars)[j], paste("beta", refcat ,sep="_")] <- unname(beta[1])
# #     results[results$inputvar==colnames(vars)[j], "p_Overall"] <- p_f
# #     results[results$inputvar==colnames(vars)[j], paste("pvalue", "_", refcat, sep="")] <- p[1]
# #     
# #     
# #     for(n in 2:ngrp){
# #       results[results$inputvar==colnames(vars)[j], paste("beta", group[n] ,sep="_")] <- unname(beta[n])
# #       results[results$inputvar==colnames(vars)[j],paste("CI_lower", refcat, sep="_")] <- unname(CIlower[1])
# #       results[results$inputvar==colnames(vars)[j],paste("CI_upper", refcat ,sep="_")] <- unname(CIupper[1])
# #       results[results$inputvar==colnames(vars)[j],paste("CI_lower", group[n] ,sep="_")] <- unname(CIlower[n])
# #       results[results$inputvar==colnames(vars)[j],paste("CI_upper", group[n] ,sep="_")] <- unname(CIupper[n])
# #       
# #       results[results$inputvar==colnames(vars)[j], paste("n",n, sep="")] <-  sum(X$group==levels(X$group)[n])
# #       results[results$inputvar==colnames(vars)[j], paste("pvalue", group[n], sep="_")] <- p[n]
# #     }
#     return(data.frame(
#       group = lm_out$xlevels$group,
#       N = sum_cases,
#       beta = dummy.coef(lm_out)$group,
#       CI = confint.lm(lm_out)[1:ngrps, ],
#       pvalue = coef(test.out)[1:ngrps,4],
#       p_overall = rep(p_f, ngrps)))
#     
#     
#   } else if (test == "fisher") {
#     Y <- matrix(c(sapply(unique(groups), function(x) sum(is.na(X[X$group==x, "y"]))), 
#                   sapply(unique(groups), function(x) sum(!is.na(X[X$group==x, "y"]))))
#                 , byrow=TRUE, ncol=ngrps, nrow=2, dimnames = list(missing = c("NA", "non-NA")))
#     ft <- try(fisher.test(Y, simulate.p.value = TRUE, B = 1e5), silent = TRUE)
#     p <- if(inherits(ft, "try-error")) {
#       rep(NA, ngrps)
#     } else {
#       ft$p.value
#     }
#     
#     sum_nas <- sapply(unique(groups), function(grp) 
#       sum(groups == grp & is.na(y)))
#     
#     
#     return(data.frame(
#       variable = unique(groups),
#       N = sum_cases,
#       NAs = sum_nas,
#       pvalue = rep(p, ngrps)))
#     
#     
#     #     results[results$inputvar==colnames(vars)[j],]$CIlower <- ft$conf.int[1]
#     #     results[results$inputvar==colnames(vars)[j],]$CIupper <- ft$conf.int[2]
#   }
# }

