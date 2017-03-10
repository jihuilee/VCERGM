#################################################
## Simulation                                   #
## Step 1) Simulate a sequence of network       #
## Step 2) Calculate the estimates of           #
##         a) Cross-sectional ERGMs             #
##         b) ad hoc 2-step procedure           #
##         c) VCERGM                            #
## Step 3) Hypothesis testing for heterogeneity #
#################################################

# Load the packge 'devtools' to install the package 'VCERGM'
library(devtools)
devtools::install_github('jihuilee/VCERGM')

library(VCERGM)

# Model with edges and mutual (reciprocity)
# Network size (nodes): n = 30 
# Number of time points: T = 50
# Network size for each observed time points: num.nodes.T = rep(n, T) 
#                                             (same for all observed networks)
# Number of simulations: nsim = 1
# Number of iterior knots for defining basis functions: interior.knot = 20
# Degrees of splines for basis functions: degree.spline = 3
# (Un)Directed networks: directed = c(TRUE, FALSE)
# Seed number when simulating networks: seed = 123

# Two network statistics (edges and reciprocity) in the model
# Use the names of statistics specified in R package 'ergm'
# Use the same strategy to specify the model 'object' as in R package 'ergm'
object = net ~ edges + mutual
K = 50
n = 30
num.nodes.K = rep(n, K) # Number of nodes at each time point
degree.spline = 3
interior.knot = 20
directed = TRUE # Directed network
nsim = 1 # Simulate one sequence of networks
seed = 123

# True phi(t) under VCERGM
phi1 = sin(((1:K) + 20)/5) + 1 # phi(t) for edges
phi2 = sin(((1:K) + 20)/3) * 0.6 + 0.4 # phi(t) for reciprocity
phi = rbind(phi1, phi2)

# True phi(t)
plot(1:K, phi1, type = "l", ylim = range(phi), xlab = "Time", ylab = "Phi", main = "True phi(t)")
lines(1:K, phi2, col = 2)

# Simulate networks
network = simulate_vcergm(object = object, num.nodes.K = num.nodes.K, phi = phi, 
                          nsim = nsim, seed = seed, directed = directed)$Networks[[1]]

# The function simulate.vcergm produces a list of simulated networks (Networks) of length nsim.
# The i-th simulated network sequence is Networks[[i]] for i = 1, ..., nsim.

# Cross-sectional ERGMs and ad hoc 2-step procedure
ergmest = cross_sectional_ergm(object = object, networks = network, directed = directed, 
                               degree.spline = degree.spline, interior.knot = interior.knot)

# cross-sectional ERGMs estimate
crossERGM_est = ergmest$phi.hat
# ad hoc 2-step procedure estimate 
adhoc_est = ergmest$phi.hat.smooth 

# VCERGM
vcergmest = estimate_vcergm(object = object, network = network, 
                            degree.spline = degree.spline, interior.knot = interior.knot, 
                            directed = directed, constant = FALSE)

# VCERGM estimate
vcergm_est = vcergmest$phi.hat

# Plotting the three estimates
plotting(ergmest, vcergmest)

# Hypothesis Test
NBoot = 1000

# Estimate under H0
est0 = estimate_vcergm(object = object, network = network, 
                       degree.spline = degree.spline, interior.knot = interior.knot, 
                       directed = directed, constant = TRUE) # Constant = TRUE (No temporal heterogeneity)

# Calculate test statistic
teststat = test_statistic(object = object, networks = network, 
                          phi0 = est0$phi.hat, phi1 = vcergmest$phi.hat, directed = directed)


# Test statistic
teststat

# Generate bootstrap samples to calculate p-value
boot = bootstrap_test(network, object = object, phi0 = est0$phi.hat, teststat = teststat,
                      degree.spline = degree.spline, interior.knot = interior.knot, directed = directed, 
                      NBoot = NBoot)

# p-value
boot$pvalue

