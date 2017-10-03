#' Simulate new networks at each time point
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'. 
#' @param num.nodes.K Number of nodes at each time point (K x 1 vector)
#' @param phi TRUE phi(t) functions (number of network statistics x K vector) 
#' @param phicoef TRUE basis coefficients (p x q matrix) 
#' @param B A set of basis functions (K x q matrix)
#' @param nsim Number of simulated networks. Default is 100.
#' @param MCMC.burnin MCMC burnin sample size. Default is 10000.
#' @param MCMC.interval Interval between selected networks. Default is 1000. The first simulated network is the (MCMC.burnin + MCMC.interval). Thereafter, every (MCMC.interval)th network will be sample.
#' @param seed Seed number used to simulate networks. Default is 123.
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#'
#' @importFrom splines bs
#' @importFrom network network
#' @importFrom network as.matrix.network
#' @importFrom ergm simulate.ergm
#' @importFrom ergm control.simulate
#' @importFrom ergm ergmMPLE
#' @export

simulate_vcergm = function(object, num.nodes.K, phi = NULL, phicoef = NULL, B = NULL,
                           nsim = 100, MCMC.burnin = 10000,
                           MCMC.interval = 1000, seed = 123, directed = c(TRUE, FALSE))
  {
  
  set.seed(seed)
  
  #  Adjust the input according to DEFAULT if not specified
  directed = directed[1]
  
  # Check that input is OK
  # 1) Check that object is a formula object
  if (class(object) != "formula" & class(object) != "ergm") {
    stop("argument object must be either a formula or ergm object. Please respecify.") }
  
  # 2) If B is NULL, replace with K x q matrix with cubic splines (assumed automatically)
  if (is.null(phi) == FALSE) {q = ncol(phi)} else{q = ncol(phicoef)} #degrees of freedom for functional basis
  K = length(num.nodes.K) # length of time series
  
  # Network statistics
  stat = unlist(strsplit(deparse(object[[3]]), " "))
  stat = stat[!stat %in% c("+", "=", "TRUE)", "FALSE)")]
  nstat = length(stat)
  
  if (is.null(B) == TRUE) {
    B = bs(1:K, df = q, degree = 3, intercept = TRUE) }

  #Initialize the list of simulated networks and the statistics for each
  network.sims = vector("list", nsim)
  
  if (is.null(phi) == FALSE) {p = nrow(phi)} else{p = nrow(phicoef)}  # number of network statistics
  
  h.statistics = rep(list(matrix(NA, K, p)), nsim)
  
  # Go through each time point and simulate nsims networks for each
  for (s in 1:K)
  {
    cat("Simulating networks for time", s, "\n")
    num.nodes = num.nodes.K[s]
    Bs = matrix(B[s,])
    
    # beta %*% B (=theta) evaluated at time point u
    if (is.null(phi) == FALSE) {coefs = as.vector(phi[,s])} else {coefs = as.vector(phicoef %*% Bs)}
    
    nets = network(num.nodes, directed = directed)
    
    #replacing object with current network formula
    z = deparse(object[[3]])
    
    formula.s = as.formula(paste("nets ~ ", z, sep = ""))
    
    # Use an existing function in package 'ergm'
    sims = simulate(object = formula.s, coef = coefs, nsim = nsim, seed = seed, 
                    control = control.simulate(MCMC.burnin = MCMC.burnin,
                                               MCMC.interval = MCMC.interval))
    
    netarray = array(NA, dim = c(num.nodes, num.nodes, nsim))
    
    if (nsim == 1)
    {
      network.sims[[1]][[s]] = as.matrix.network(sims, matrix.type = "adjacency")
      h.statistics[[1]][s,] = summary(as.formula(paste("sims ~ ", z, sep = "")))
#      if (is.null(rownames(phicoef)) == FALSE) {colnames(h.statistics[[1]]) = rownames(phicoef)}
      colnames(h.statistics[[1]]) = stat
    }
    
    if (nsim > 1)
    {
      for (i in 1:nsim)
      {
        network.sims[[i]][[s]] = as.matrix(sims[[i]], matrix.type = "adjacency")
        h.statistics[[i]][s, ] = as.matrix(attr(sims, "stats"))[i, ]
#        if (is.null(rownames(phicoef)) == FALSE) {colnames(h.statistics[[i]]) = rownames(phicoef)}
        colnames(h.statistics[[i]]) = stat
      }
    }
    
  }
  return(list(Networks = network.sims, Statistics = h.statistics))
}

