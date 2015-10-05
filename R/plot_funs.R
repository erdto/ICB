

#' Plotting function for pipeline_pairwise
#' 
#' takes data.frame with results (from compare_grps function)
#' takes data.frame without NAs
#' @export

print_boxplot <- function(X, pvalue, name, boxplot_args) {
  stopifnot(is.data.frame(X) && "group" %in% names(X))
  stopifnot(all(names(X) %in% c("y", "group")))
  ggplot <- boxplot_args$ggplot
  
  X <- na.exclude(X)
  title <- sprintf("p = %.6g", pvalue)
  
  # if there are no covariates, var has NAs and we have to set those
  # same rows NA in the grp variable
  #   grp <- group
  #   grp[which(is.na(var))] <- NA
  
  # when there are covariates, var is already cleansed from NAs
  # when residuals were computed. So we have to only take rows from grp
  # that correspond to the row names: names(var)
  #   if(length(grp) != length(var)) {
  #     grp <- grp[as.numeric(names(var))]
  #   }
  
  if(ggplot) {
    theme <- ggplot2::theme_set(ggplot2::theme_minimal())
    q <- ggplot2::ggplot(data = X) +
    ggplot2::geom_boxplot(ggplot2::aes(x = group, y = y, colour = factor(group))) +
      ggplot2::ylab(name) +
      ggplot2::geom_jitter(ggplot2::aes(x = group, y = y, color = factor(group))) +
      ggplot2::ggtitle(title) +
      ggplot2::scale_x_discrete("Group", labels = levels(X$group)) + #old  c(compare[1], compare[2])
      ggplot2::theme(legend.position="none") 
    
    if(!is.null(boxplot_args$ylim))
      q <- q + ggplot2::ylim(boxplot_args$ylim)
    
    print(q)
    return(q)
  } else {
    boxplot(X[,1] ~ group, data = X, col = "white", border = c("coral2", "turquoise3"), lwd = 3,
            main = title, ylab = name)
    stripchart(X[,1] ~ group, vertical = TRUE, data = X, 
               method = "jitter", add = TRUE, pch = 16, col = c("coral2", "turquoise3"))
  }
}

###########################################################################################
##############

#' Plotting function for pipeline_overall
#'  @export

print_boxplot2 <- function(X, name, p_overall, boxplot_args) {
  stopifnot(is.data.frame(X) && all(c("y", "group") %in% names(X)))
  stopifnot(is.list(boxplot_args))
  
  if(all(is.na(X[,1])))
    return(NULL)
  
  ggplot <- boxplot_args$ggplot
  boxplot_colors = boxplot_args$boxplot_colors
  boxplot_fill = boxplot_args$boxplot_fill
  boxplot_size = boxplot_args$boxplot_size
  boxplot_shape = boxplot_args$boxplot_shape
  boxplot_alpha = boxplot_args$boxplot_alpha
  
  # get colors
#   if(length(unique(group)) != boxplot_colors) {
#     hues <- seq(15, 375, length=number_vals+1)
#     cols <- hcl(h=hues, l=65, c=100)[1:number_vals]
#   }
  
  title <- sprintf("Overall test p = %.6g", p_overall)
  X <- na.exclude(X) 
  
  if(ncol(X) > 2)
    X$yb <- lm(X[ , "y"] ~ . -group, data=X)$res
  else
    X$yb <- X[ , "y"]
  
  group <- X$group 
  if(ggplot) {
    theme <- ggplot2::theme_set(ggplot2::theme_minimal())
    b <- ggplot2::qplot(factor(group), yb, data=X, geom="boxplot", colour=factor(group)) + 
      ggplot2::ylab(name) +  ggplot2::ggtitle(title) +
      ggplot2::geom_jitter(ggplot2::aes(color = factor(group))) + 
      ggplot2::scale_x_discrete("Group", labels = levels(group)) + 
      ggplot2::theme(legend.position="none") + 
      # ggplot2::scale_colour_manual(values=boxplot_colors) +
      if(!is.null(boxplot_args$ylim))
        ggplot2::ylim(ylim) +
      ggplot2::scale_size_manual(values= boxplot_size)  +
      ggplot2::scale_alpha_manual(values = boxplot_alpha) +
      ggplot2::scale_shape_manual(values = boxplot_shape) +
      ggplot2::scale_fill_manual(values = boxplot_fill)
    print(b)
  } else {
    ## without ggplot
    boxplot(yb~group,data=X, main=title, lwd = 2,
            xlab="Group", ylab=names(X)[1],  col = boxplot_colors)
  }
}

###########################################################################################
##############

#' Takes data.frame with results from compare_grps function
#' and volc_args: list with parameters that is set by function set_volc_args()
#' prints and returns the plot
#'  @export
print_volc <- function(datframe, volc_args) {
  stopifnot(is.data.frame(datframe) & all(c("beta", "pvalue", "inputvar") %in% names(datframe)))
  stopifnot(all(
    c("volc_cutoff", "volc_colors", "volc_size", "volc_shape", "volc_alpha") 
    %in% names(volc_args)))
  
  volc_cutoff = volc_args$volc_cutoff
  volc_colors = volc_args$volc_colors
  volc_size = volc_args$volc_size
  volc_shape = volc_args$volc_shape
  volc_alpha = volc_args$volc_alpha
  
  # set ylim from whole data
  ylim <- c(0, max(-log10(abs(datframe$pvalue)), na.rm=T)+5)
  xlim <- c(-max(datframe$beta, na.rm=T)-5, max(datframe$beta, na.rm=T)+5)
  
  # use only values with a non-NA beta coefficient
  datframe <- datframe[with(datframe, !is.na(beta)),]
  
  # Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
  datframe$threshold <- as.factor(abs(datframe$beta) > volc_cutoff["beta"] &
                                    abs(log10(datframe$pvalue)) > -log10(volc_cutoff["p"]))
  
  # delete labels of sub-threshold objects
  cols <- c("FALSE" = volc_colors[1], "TRUE" = volc_colors[2])
  # Construct the plot object
  ggplot2::theme_set(ggplot2::theme_minimal())
  suppressWarnings(
    g <- ggplot2::ggplot(data = datframe, ggplot2::aes(x = beta, y= -log10(pvalue),
                                                       colour = threshold, fill = threshold)) +
      ggplot2::geom_point(alpha = volc_alpha, size = volc_size, shape = volc_shape) +
      ggplot2::xlim(xlim) +
      ggplot2::ylim(ylim) +
      ggplot2::xlab(paste("<- log fold change ->")) +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::scale_colour_manual(values = cols) +
      ggplot2::scale_fill_manual(values = cols)
  )
  # annotate points beyond the significance lines
  set_na <- !is.na(datframe$pvalue) & datframe$threshold == TRUE
  
  #     union(which(datframe$threshold == FALSE), which(is.na(datframe$pvalue)))
  datframe[!set_na, "inputvar"] <- NA
  if(any(!is.na(datframe$inputvar))) {
    datframe[datframe$beta < 0, "hj"] <- 1
    datframe[datframe$beta > 0, "hj"] <- 0
    suppressWarnings(g <- 
                       g + ggplot2::geom_text(data = datframe, 
                                              ggplot2::aes(x = beta, 
                                                           y = -log10(pvalue),
                                                           label = inputvar, size = 1.0,  hjust = hj), 
                                              vjust = 0, colour = "black")
    )
  }
  
  # add 0 and "significance" lines
  suppressWarnings(g <- 
                     g  +  ggplot2::theme(legend.position="none") +
                     ggplot2::geom_vline(xintercept=0, linetype="longdash", colour="black") +
                     ggplot2::geom_hline(yintercept=(-log10(volc_cutoff["p"])), linetype="dotted", colour="black") +
                     ggplot2::geom_vline(xintercept = volc_cutoff["beta"],  linetype="dotted", colour="black") +
                     ggplot2::geom_vline(xintercept=(-volc_cutoff["beta"]), linetype="dotted", colour="black")
  )
  
  #   if(!is.null(volc_args)){g <- g + volc_args}
  suppressWarnings( print(g) )
  g
}

###############################################################################
#########

#' @export
print_scatterplot <- function(y, predictor, covariates, name, predname) {
  
  # partial-out effects of covars; compute and save partial correlation and y-residuals for plotting
  if(nrow(covariates) > 0) {
    df <- data.frame(y = y, x = predictor, covariates)
    df <- na.exclude(df)
    yres <- lm(y~ ., data = df[, -2])$res
    xres <- lm(x~., data = df[, -1])$res
  } else {
    df <- data.frame(y = y, x = predictor)
    df <- na.exclude(df)
    yres <- df$y
    xres <- df$x
  }
  
  X <- data.frame(yres = yres, xres = xres)
  if(nrow(X) < 3)
    return(NULL)
  
  partialcor <- cor(yres, xres)
  
  title <- sprintf("pcor = %.6g", partialcor)
  theme <- ggplot2::theme_set(ggplot2::theme_minimal())
  
  b <- ggplot2::qplot(xres, yres, data=X, geom="point", colour=xres) + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab(name) + 
    ggplot2::xlab(predname) + 
    ggplot2::ggtitle(title) +
    ggplot2::theme(legend.position="none")
  
  print(b)
}

###############################################################################
#########


#' Takes metabolite data and covariates
make_boxplots <- function(results, ...) UseMethod("make_boxplots")


make_boxplots.pairwise <- function(results, vars, groups, covariates, boxplot_args) {
  stopifnot(all(c("sort_boxplots", "ylim", "ggplot") %in% names(boxplot_args)))
  
  # sort boxplots after pvalues
  if(boxplot_args$sort_boxplots) {
    neworder <- order(results$pvalue, decreasing = FALSE)
    results <- results[neworder, ]
    vars <- vars[, neworder]
  } 
  
  checkmate::assertDataFrame(x = vars, all.missing = FALSE, types = "numeric", min.cols = 1, min.rows = 2)
  # compute residuals, i.e. partial out effects of covariates
  if(NROW(covariates) == 0) {
    residuals <- lapply(vars, function(v) {
      names(v) <- 1:length(v)
      return(v) 
      })
  } else {
    
    # get residuals to order to partial out the effects of the covariables before plotting the 
    # group effects.
    # If this fails, it will be because of NAs leading to some factors having only one unique value,
    # in this case, just print the unpartialed effects.
    residuals <- setNames(sapply(1:ncol(vars), function(ii) {
      y <- vars[, ii]
      if(is.null(names(y)))
        y <- setNames(y, nm = 1:length(y))
      rval <- try(lm(y ~., data = data.frame(cbind(y, covariates)))$res, silent =  TRUE)
      if(inherits(rval, "try-error"))
        return(vars[,ii])
      rval
    }, simplify = FALSE), nm = names(vars))
    
  }
  
  plotlist <- lapply(names(residuals), function(x) {
    stopifnot(all(as.numeric(names(residuals[[x]])) %in% 1:length(groups)))
    if(!is.na(results[results$inputvar == x, "pvalue"]))
      print_boxplot(X = na.exclude(data.frame(y = residuals[[x]], 
                                              group = groups[as.numeric(names(residuals[[x]]))])),
                    name = x,
                    pvalue = results[results$inputvar == x, "pvalue"],
                    boxplot_args = boxplot_args) 
  })
  plotlist
}

make_boxplots.overall <- function(results, vars, groups) {
  checkmate::assertDataFrame(x = vars, all.missing = FALSE, types = "numeric", min.cols = 1, min.rows = 2)
  
  # TODO, partial out covars!
  
  # boxplot_args
  sapply(names(vars), function(X) {
    if(!is.na(results[results$var == X, "p_overall"][1]))
      print_boxplots2(data.frame(X = vars[, X],
                                 group = groups),
                      p_overall =  results[results$var == X, "p_overall"][1], 
                      boxplot_args = set_boxplot_args2())})
  
}

make_boxplots.continuous <- function(results, vars, predictor, covariates) {
  
  checkmate::assertDataFrame(x = vars, all.missing = FALSE, types = "numeric", 
                             min.cols = 1, min.rows = 2)
  
  
  if(nrow(covariates) > 0) {
    X <- data.frame(y = )
    lm1 <- lm(y ~ X[, 3:ncol(X)], data=X)
    lm2 <- lm(pred ~ X[, 3:ncol(X)], data=X)
    
  }
  #               # partial-out effects of covars; compute and save partial correlation and y-residuals for plotting
  #               X <- na.exclude(X) 
  #               
  #                 lm1 <- lm(y ~ X[, 3:ncol(X)], data=X)
  #                 lm2 <- lm(X[,2] ~ X[, 3:ncol(X)], data=X)
  #               X$yres <- lm1$res
  #               X$xres <- lm2$res
  #   partialcor <- sapply(vars, function(j))
  #                 results[results$inputvar==colnames(vars)[j], "partialcor"] <- cor(X$yres, X$xres, method = cor_method)
  #               
  #               title <- sprintf("pcor = %.6g", results[results$inputvar==colnames(vars)[j], "partialcor"])
  #               theme <- theme_set(theme_minimal())
  #               b <- qplot(X[,2], yres, data=X, geom="point", colour=X[,2]) + geom_smooth(method = "lm") +
  #                                 ylab(colnames(vars)[j]) + xlab(names(X)[2]) + ggtitle(title) +
  #                                 theme(legend.position="none") #+ 
  #               #                 scale_colour_manual(values=boxplot.colours) +
  #               #                 scale_size_manual(values= boxplot.size)  +
  #               #                 scale_alpha_manual(values = boxplot.alpha) +
  #               #                 scale_shape_manual(values = boxplot.shape) +
  #               #                 scale_fill_manual(values = boxplot.fill)
  #               
  #               if(class(boxplot_theme)[1]=="theme"){b <- b + boxplot_theme}  #+ boxplot_layers #theme_update(theme, unlist(boxplot_theme))
  #               
  #               print(b)
  #             }
  #           }
  #         }) # vars loop over
  
  #       R <- results[!is.na(results$beta), ]  
  #       R <- R[is.finite(-log10(results$pvalue)), ]
  # this was for closing scatterplot file
  #           dev.off()
  
}





