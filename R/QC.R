#' Allele Specific Expression Quality Control (ASEQC) Test
#'
#' This test is designed to detect poor quality ASE data. This includes
#' low coverage samples, samples with potential genotype mismatches while
#' generating the ASE data etc.
#'
#' @include blnmix.R
#' @include preprocess.R

#' @import parallel
#' @import doParallel
#' @import foreach
#' @import pracma
#' @import dplyr
#' @import stringr
#' @import data.table
#' @import R.utils
#' @import robustbase
#' @import stats
#' @import utils
#'
#' @param ref Dataframe containing reference count data
#' @param alt Dataframe containing alternate count data
#' @param numCores number of cores used for parallel execution (default = 1).
#'
#'
#' @return  Fitted std dev. and the ascertained quality cutoff.
#' @export


QC <- function(ref, alt, numCores = 1) {
  # input validation
  if (!is.numeric(numCores)) {
    stop("Please enter a number for number of cores")
  }
  if (!is.data.frame(ref) || !is.data.frame(alt)) {
    stop("Either ref or alt is not of type dataframe")
  }
  if (ncol(ref) != ncol(alt) || nrow(ref) != nrow(alt)) {
    stop("Datasets contain different dimensions")
  }
  # Initializing dataframe
  samples_output_list <- colnames(ref)[-2]

  # Temporary Files
  file_out <- tempfile(pattern = "file.out.txt")
  file_indices <- tempfile(pattern = "file.indices.txt")

  # Remove Temporary Files if existant.
  if (file.exists(file_out)) {
    file.remove(file_out)
  }

  if (file.exists(file_indices)) {
    file.remove(file_indices)
  }

  # For the parallel loop
  cluster_for_loop <- makeCluster(numCores)
  registerDoParallel(cluster_for_loop)

  # Foreach loop doing things parallely
  foreach(i = 1:ncol(ref), .packages = "pracma", .combine = rbind) %dopar% {

    # Getting Sample
    if (!is.numeric(ref[, i]) || !is.numeric(alt[, i])) {
      writeLines("Non numerical Row: Moving to next Row")
    } else {

      # Here, r refers to the reference allele counts.
      # And a, refers to the alternate allele counts.
      # t, refers to the total allele counts, sum of r and a

      r <- (ref[, i])
      a <- (alt[, i])
      t <- r + a

      # First remove the total counts which do not meet the ascertained
      # total_count threshold.
      # Here, we return those total counts that are between 5 and 5000

      lst <- remove_total(r, t)

      # Fit the modified BLN distribution to the data
      params <- Fit_BLN_uniform_mixture(lst$ref, lst$total)
      # Saving to vectors
      params <- c(params, lst$num.kept, median(lst$total))

      # Write to file, this is to save the order in which
      # the input file was written.

      # When  each sample is sent to a number of cores
      # The second line of code is to keep track which sample
      # index in the original file it was.

      write(params, file_out, append = TRUE)
      write(i, file_indices, append = TRUE)
    }
  }
  # Once done, stop the cluster.
  stopCluster(cluster_for_loop)

  # Reading written to recreate the original table.
  df_output <- read.table(file_out, sep = " ", fill=TRUE)
  sample_indices <- read.table(file_indices, sep = " ", fill=TRUE)

  # Setting the file index
  df_output$indices <- sample_indices$V1

  # Sorting in Order of original input
  df_output <- setorder(df_output, indices)

  # Adding in Samples
  df_output$sample <- samples_output_list

  # Changing column names
  colnames(df_output) <- c("mean", "std", "num.kept", "median.coverage", "indices", "sample")

  cols_to_keep <- c("mean", "std", "num.kept", "median.coverage", "sample")

  df_output <- df_output[, cols_to_keep]

  # Remove fit for index row
  df_output <- df_output %>% slice(-1)

  # Add in the skew-adjusted boxplot upper whisker.
  threshold <- adjboxStats(df_output$std)$stats[5]
  df_output$cutoff <- threshold

  return(df_output)
}
