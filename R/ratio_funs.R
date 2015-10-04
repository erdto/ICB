#----------------------------------------------------------------------------------
### Utility functions for pipelines

#----------------------------------------------------------------------------------
# takes data.frame and column with grouping var
# and returns data.frame with new columns with aggregated z scores
use_agg_z_df <- function(df, vargrp){
  
  varcomb <- unique(vargrp)[!is.na(unique(vargrp))]
  #     newdata <- data.frame(matrix(NA, nrow = nrow(df), ncol = length(varcomb)))
  #     names(newdata) <- varcomb
  newdat <- as.data.frame(lapply(varcomb, function(x) {
    new_columns <- which(vargrp %in% x)
    apply(scale(df[ , new_columns[which(new_columns <= ncol(df))]]),
          1, function(x) mean(x, na.rm=T))
  }))
  names(newdat) <- varcomb
  newdat
}

###########################################################################################
##############

#' df data.frame with columns corresponding to metabolites (or similar).
#' names [vector] of metabolite names.
#' ratioop one of c("subtract", "divide") -> use subtract to make ratios for logarithmic data.
add_ratios_to_df <- function(df, names, ratioop = c("subtract", "divide")) {
  ratioop <- match.arg(ratioop)
  ratio <- combn(names,2)
  
  if(ratioop=="divide") {
    vratios <- sapply(1:ncol(ratio),  function(x)
      df[ ,which(names %in% ratio[1, x])] / df[ , which(names %in% ratio[2, x])], simplify = FALSE )
  } else if (ratioop == "subtract") {
    vratios <- sapply(1:ncol(ratio),  function(x)
      df[ ,which(names %in% ratio[1, x])] - df[ , which(names %in% ratio[2, x])], simplify = FALSE )
  }
  names(vratios) <- paste(ratio[1,], "/",ratio[2,], sep="")
  cbind(df, vratios)
}

###########################################################################################
##############

# compute_pgain2 <- function(R) {
#   stopifnot(is.data.frame(R) && all(c("pvalue", "inputvar", "denominator") %in% names(R)))
#   
#   
#   name1 <- R$inputvar
#   name2 <- R$denominator
#   ratio_inds <- which(!is.na(R$denominator))
#   
#   pratio <- R[ratio_inds, "pvalue"]
#   p1 <- unlist(sapply(name1[ratio_inds], function(x)
#     R[R$inputvar == x & is.na(R$denominator), "pvalue"]))
#   p2 <- unlist(sapply(name2[ratio_inds], function(x)
#     R[R$inputvar == x & is.na(R$denominator), "pvalue"]))
#   
#   pgain <- pratio / pmin(p1, p2)
#   rval <- numeric(nrow(R))
#   rval[!is.na(R$pvalue)] <- pgain
#   rval
# }

#' add p-gain column to output data.frame
#' R is data.frame with results of the tests.
#' Returns data.frame with additional column.
#' 
compute_pgain <- function(R) {
  stopifnot(is.data.frame(R) && all(c("pvalue", "inputvar", "denominator") %in% names(R)))
  
  pgain <- numeric(nrow(R))
  
  pvals <- R[which(is.na(R[, "denominator"])), c("inputvar", "pvalue")]
  ratioinds <- which(!is.na(R[ , "denominator"]))
  
  for(i in ratioinds) {
    pratio <- R[i, "pvalue"]
    p1 <- pvals[pvals$inputvar == R[i, "inputvar"], "pvalue"]
    p2 <- pvals[pvals$inputvar == R[i, "denominator"], "pvalue"]
    pgain[i] <- pratio / min(p1, p2)
  }
  pgain
}


###########################################################################################
##############
