#' Extracts total counts from reference and alternative counts
#'
#' @param ref Reference read counts
#' @param alt Alternative read counts
#' @return total counts as a dataframe
#' @noRd

extract_total <- function(ref, alt) {
  total <- ref[,2:ncol(ref)] + alt[,2:ncol(alt)]
  total <- cbind('name'=ref$name, total)
  return(total)
}

#' Removes Expression less or greater than the given thresholds to speed
#' up optimization.
#'
#' @param ref Reference read counts.
#' @param total Total read counts.
#' @param min Threshold to remove total read counts less than
#' @param max Threshold to remove total read counts greater than
#' @return Updated reference counts, total counts, and number of counts removed
#' @noRd

remove_total <- function(ref, total, min=5, max=5000) {
  vec <- c()
  for(i in seq(length(total))) {
    if(((total[i] <= min) || (total[i] >= max))) {
      vec <- c(vec,-i)
    }
  }
  if(length(vec) == 0) { # Case where we removed nothing
    lst <- list(ref = ref, total = total,removed = 0)
    return(lst)
  }
  lst <- list(ref = ref[vec], total = total[vec], num.kept = length(total[vec]))
  return(lst)
}
