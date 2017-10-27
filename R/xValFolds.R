#' xValFolds function taken from cvTools
#'
#' @note This function is taken from the cvTools package.
#'
#' @concept rcellminerElasticNet
#' @export
#' 
xValFolds <- function (n, K = 5, R = 1, type = c("random", "consecutive", 
                                    "interleaved")) 
{
  n <- round(rep(n, length.out = 1))
  if (!isTRUE(n > 0)) 
    stop("'n' must be positive")
  K <- round(rep(K, length.out = 1))
  if (!isTRUE((K > 1) && K <= n)) 
    stop("'K' outside allowable range")
  type <- if (K == n) 
    "leave-one-out"
  else match.arg(type)
  if (type == "random") {
    R <- round(rep(R, length.out = 1))
    if (!isTRUE(R > 0)) 
      R <- 1
    subsets <- replicate(R, sample(n))
  }
  else {
    R <- 1
    subsets <- as.matrix(seq_len(n))
  }
  which <- rep(seq_len(K), length.out = n)
  if (type == "consecutive") 
    which <- rep.int(seq_len(K), tabulate(which))
  folds <- list(n = n, K = K, R = R, subsets = subsets, which = which)
  #class(folds) <- "cvFolds"
  folds
}