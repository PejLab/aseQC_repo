#' Fit_BLN_uniform_mixture
#'
#' @param ref The reference counts for the sample by gene
#' @param total The total counts for the sample by gene
#' @param start An array containing the start points for the optimisation.
#' @importFrom bln dbln
#' @return The mean and std dev. of a BLN-Uniform mixture distribution.
#' @noRd

Fit_BLN_uniform_mixture <- function(ref, total, start) {

  if(missing(start)) {
    start.Std <- 1 # starting point for sigma prediction
    start.Mean <- 0 # Starting point for mu prediction

  } else {

    start.Mean <- start[1]
    start.Std <- start[2]
  }

  MLE_optim <- optim(par <- c(start.Mean, start.Std), fn = BLN_uniform_mixture_likelihood2,
                     ref = ref, total = total, lower = c(-0.25,0), upper = c(0.25,Inf),
                     method = "L-BFGS-B")


  return(MLE_optim$par)
}


#' BLN_uniform_mixture_likelihood2
#'
#' @param params An array of mean and std dev for the BLN.
#' @param ref The reference counts for the sample by gene
#' @param total The total counts for the sample by gene
#' @importFrom bln dbln
#'
#' @return NLL, the negative log likelihood of the BLN with input parameters
#'
#' @noRd
#'
BLN_uniform_mixture_likelihood2 <- function(params,ref, total) {

  # parameters
  mean <- params[1]
  sd <- params[2]

  # Calculating Likelihood
  BLN_L <- dbln(ref, total, mean, sd)
  UNIF_L <- (1/(total + 1))
  LIKE <- .999 * BLN_L + (1 - .999) * UNIF_L
  NLL <- -sum(log(LIKE))
  return(NLL)
}
