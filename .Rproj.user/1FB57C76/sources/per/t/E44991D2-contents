#' Differentially Expressed Heterogeneous Overdispersion Genes Testing for Count Data
#' This script implements the main function of the proposed method in the above paper

#' @param data A matrix of gene expression data where rows represent genes and columns represent samples.
#' @param treatment A vector specifying the treatment conditions for each sample.
#' @param covariates An optional matrix of gene-wise covariates. Default is NULL.
#' @param norm_factors An optional vector of normalization factors for each sample. Default is NULL, which assumes equal normalization factors.
#' @param dist The distribution family for the GLM. Can be "qpois" for quasi-Poisson or "negbin" for negative binomial. Default is "qpois".
#' @param padj Logical value indicating whether to adjust p-values using the Benjamini-Hochberg (BH) procedure. Default is TRUE.
#' @param pval_thre The threshold for identifying differentially expressed genes based on adjusted p-values. Default is 0.05.
#' @param l2fc Logical value indicating whether to consider log2 fold change for identifying differentially expressed genes. Default is FALSE.
#' @param l2fc_thre The threshold for log2 fold change in identifying differentially expressed genes. Default is 1.
#' @param num_cores The number of CPU cores to use for parallel computing. Default is 1.
#' @return A list containing:
#'   \item{DE_idx}{A logical vector indicating differentially expressed genes.}
#'   \item{pvals}{A numeric vector of p-values for each gene.}
#'   \item{log2fc}{A numeric vector of log2 fold changes for each gene.}
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom MASS glm.nb
#' @importFrom stats as.formula coef glm median p.adjust quasipoisson
#' @examples
#' # simulate gene expression data
#' data <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' # simulate random treatment assignments
#' treatment <- sample(0:1, 10, replace = TRUE)
#' # Run main function with parallel computing using 2 cores
#' result <- dehogt_func(data, treatment, num_cores = 2)
#' @export
dehogt_func <- function(data, treatment, norm_factors = NULL, covariates=NULL, dist = "qpois", padj = TRUE, pval_thre = 0.05, l2fc = FALSE, l2fc_thre = 1, num_cores = 1) {

  median_ratio_norm_factor <- function(data){

    # create a pseudo-reference sample (geometric mean across all samples for each gene)
    geom_means <- apply(data, 1, function(x) exp(mean(log(x[x > 0]))))

    # calculate the ratio of each sample to the pseudo-reference sample
    ratios <- t(t(data) / geom_means)

    # Step 3: Calculate the normalization factor (median of ratios for each sample)
    size_factors <- apply(ratios, 2, median, na.rm = TRUE)

    return(size_factors)
  }

  # Initialize parallel computing
  registerDoParallel(cores = num_cores)

  # number of samples
  num_samples <- ncol(data)

  # covariates adjustment if covariates are provided

  if (!is.null(covariates)) {

    if (nrow(covariates) != nrow(data)) {
      stop("Number of rows in covariates matrix should match the number of genes in the data matrix.")
    }

    else {

      vec_data = as.vector(data)

      rep_row_index = rep(1: nrow(data), each = num_samples)

      rep_covariates = covariates[rep_row_index, ]

      if (dist == "qpois") {
        adjust_model = glm(vec_data ~ rep_covariates - 1, family = quasipoisson())
      } else if (dist == "negbin") {
        adjust_model = glm.nb(vec_data ~ rep_covariates - 1)
      }

      beta = coef(adjust_model)
    }
  }

  # calculate size factors for each sample
  size_factors <- median_ratio_norm_factor(data)

  # perform gene-wise GLM, with parallel computing
  result_GLM <- foreach(i = seq_len(nrow(data)), .combine = rbind) %dopar% {

    genewise_data <- data[i, ]

    if (!is.null(covariates)) {
      adjust <- covariates[i, ] %*% beta
    } else {
      adjust <- 0
    }

    genewise_dataframe <- data.frame(count = genewise_data, treatment = treatment, size_factors = size_factors, adjust = rep(adjust, num_samples))

    # define formula for GLM
    model_formula = as.formula("count ~ treatment + offset(log(size_factors)) + offset(adjust)")

    if (dist == "qpois") {
      model = glm(formula = model_formula, family = quasipoisson(), data = genewise_dataframe)

      coefficients <- coef(summary(model))[2, "Estimate"]
      p_values <- coef(summary(model))[2, "Pr(>|t|)"]
      padj_values <- p.adjust(p_values, method = "BH")
      log2fold <- log(exp(abs(coefficients)), base = 2)

    } else if (dist == "negbin") {
      model = glm.nb(formula = model_formula, data = genewise_dataframe)

      coefficients <- coef(summary(model))[2, "Estimate"]
      p_values <- coef(summary(model))[2, "Pr(>|z|)"]
      padj_values <- p.adjust(p_values, method = "BH")
      log2fold <- log(exp(abs(coefficients)), base = 2)

    }

    c(coefficients, p_values, padj_values, log2fold)
  }

  stopImplicitCluster()

  colnames(result_GLM) <- c("Coefficient", "P_value", "P_adj_value", "Log2fold")

  if (padj) {
    pvals <- result_GLM[, "P_adj_value"]
  } else {
    pvals <- result_GLM[, "P_value"]
  }

  idx <- which(pvals < pval_thre)

  log2fold <- result_GLM[, "Log2fold"]

  if (l2fc) {
    idx <- which(abs(log2fold) > l2fc_thre)
  }

  output <- list(DE_idx = as.numeric(idx), pvals = as.numeric(pvals), log2fc = as.numeric(log2fold))

  return(output)
}
