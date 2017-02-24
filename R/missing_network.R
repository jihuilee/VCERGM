#' Generate missing networks
#'
#' @param networks A list of (fully) observed networks
#' @param num.missing Number of missing networks. Default is NULL.
#' @param missing.time A list of time points with missing networks. Default is NULL.
#' @param seed Seed number to generate missing.time when only num.missing is provided. Default seed number is 1.

missing.network = function(networks, num.missing = NULL, missing.time = NULL, seed = 1)
  
{
  set.seed(seed)
  K = length(networks)
  
  # Make sure there is no duplicated time points.
  missing.time = unique(missing.time)
    
  ## If both num.missing and missing.time are specified 
  ## 1) with different lengths, ERROR!
  ## 2) with the same length, use missing.time.
  
  if (is.null(num.missing) != TRUE & is.null(missing.time) != TRUE)
  {
    if (length(missing.time) != num.missing)
    {stop("Both num.missing and missing.time are specified and they have different lengths. Please respecify.")}    
  }
  
  if (is.null(num.missing) != TRUE & is.null(missing.time))
  {missing.time = sample(1:K, num.missing, replace = FALSE)}
  
  if (is.null(missing.time) != TRUE & max(missing.time) > K)
  {
    stop("You include at least one time point beyond the actual length of observed network.
         Please respecify")
  }
    
  missing.time = sort(missing.time)
  
  newnetworks = list()
  for (i in 1:K)
  {
    newnetworks[[i]] = networks[[i]]
    if (i %in% missing.time) {newnetworks[[i]] = NULL}
  }
  return(list(newnetworks = newnetworks, missing.time = missing.time))
}
