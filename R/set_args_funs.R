
#' Set arguments for ratio option in \link{\code{pipeline_pairwise}}, \link{\code{pipeline_continous}} and \link{\code{pipeline_overall}}.
#' @param use_ratios [logcial] should ratios of the metabolite columns be computed and analyzed?
#' @param ratio_op [choice] are ratios computed by subtracting (for logarithmic data) or dividing?
#' @return list with arguments.
set_ratio_args <- function(use_ratios = FALSE, 
                           ratio_op = c("subtract", "divide")) {
  ratio_op <- match.arg(ratio_op)
  list(ratios = use_ratios, ratio_op = ratio_op)
}


#' Set arguments for volcano plots \link{\code{pipeline_pairwise}} and \link{\code{pipeline_continuous}}.
#' @param cutoff [numeric vector] Set cut-off values for the coefficient and p-value. Metabolites with values beyond the cut-off will be colored and labeled.
#' @param cutoff_fdr [logical] Should the FDR corrected p-values (instead of the naive, uncorrected) be used as cutoff?
#' @param colors colors
#' @param size size of the points
#' @param shape shape 
#' @param fill fillcolor
#' @param alpha transparency
#' @export
set_volc_args <- function(cutoff = c("beta" = 0.5, "p" = 0.05), 
                          cutoff_fdr = FALSE,
                          colors = c("coral2", "turquoise3"),
                          size = 1.75,
                          shape = 21,
                          fill = "white",
                          alpha = 1) {
  list(volc_cutoff = cutoff, volc_cutoff_fdr = cutoff_fdr, volc_colors = colors, volc_size = size,
       volc_shape = shape, volc_fill = fill, volc_alpha = alpha)
}

#' Set arguments for boxplots in \link{\code{pipeline_pairwise}}
#' @param sort_boxplots should boxplots be sorted after p-values?
#' @param ylim [numeric vector(length = 2)] set ylim for all boxplots.
#' @param p_cutoff [numeric] Only metabolites with below cut-off p-value will be shown in boxplots. Default is no cut-off.
#' @param ggplot should ggplot2 package be used for plotting?
#' @param colors colors
#' @param size size of the points
#' @param shape shape 
#' @param fill fillcolor
#' @param alpha transparency
#' @export
set_boxplot_args <- function(sort_boxplots = FALSE,
                             ylim = NULL,
                             p_cutoff = Inf,
                             ggplot = TRUE, 
                             colors = c("coral2", "turquoise3"), 
                             fill = "white",
                             size = 1, 
                             shape = 21, 
                             alpha = 1) {
  list(sort_boxplots = sort_boxplots, ylim = ylim, p_cutoff = p_cutoff, ggplot = ggplot, boxplot_colors = colors, boxplot_fill = fill, boxplot_size = size,
       boxplot_shape = shape, boxplot_alpha = alpha)
}

#' Set arguments for boxplots in \link{\code{pipeline_overall}}
#' @param sort_boxplots should boxplots be sorted after p-values?
#' @param no_plots [logical] Should all plots be omitted? This will speed up the function considerably.
#' @param ylim [numeric vector(length = 2)] set ylim for all boxplots.
#' @param p_cutoff [numeric] Only metabolites with below cut-off p-value will be shown in boxplots. Default is no cut-off.
#' @param colors colors
#' @param ggplot should ggplot2 package be used for plotting?
#' @param size size of the points
#' @param shape shape 
#' @param fill fillcolor
#' @param alpha transparency
#' @export
set_boxplot_args2 <- function(sort_boxplots = FALSE,
                              no_plots = FALSE,
                              ylim = NULL,
                              p_cutoff = Inf,
                              colors = c("#E69F00", "#56B4E9", "#009E73"),
                              ggplot = TRUE, 
                              fill = "white",
                              size = 1, 
                              shape = 21, 
                              alpha = 1) {
  list(sort_boxplots = sort_boxplots, no_plots = no_plots, ylim = ylim, p_cutoff = p_cutoff,
       ggplot = ggplot, boxplot_colors = colors, boxplot_fill = fill, boxplot_size = size,
       boxplot_shape = shape, boxplot_alpha = alpha)
}
