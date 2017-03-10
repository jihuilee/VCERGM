#################################################
## Application (Voting data)                    #
## Step 1) Calculate the estimates of           #
##         a) Cross-sectional ERGMs             #
##         b) ad hoc 2-step procedure           #
##         c) VCERGM                            #
## Step 2) Hypothesis testing for heterogeneity #
#################################################

# Load the packge 'devtools' to install the package 'VCERGM'
library(devtools)
devtools::install_github('jihuilee/VCERGM')

library(VCERGM)

# Load the data (Congress 40 - 113)
data(Rollcall)

networks = Rollcall$networks # Networks
attr = Rollcall$attr # Political affiliation

object = Net ~ triangle + kstar(2) + nodematch("attr1")

degree.spline = 3
interior.knot = 30
directed = FALSE

# Cross-sectional ERGMs and ad hoc 2-step procedure
ergmest = cross_sectional_ergm(object = object, network = networks, attr = attr, directed = directed, 
                     degree.spline = degree.spline, interior.knot = interior.knot)

crossERGM_est = ergmest$phi.hat
adhoc_est = ergmest$phi.hat.smooth

# VCERGM
vcergmest = estimate_vcergm(object = object, network = networks, attr = attr,
                            degree.spline = degree.spline, interior.knot = interior.knot, 
                            directed = directed, constant = FALSE)

vcergm_est = vcergmest$phicoef.hat

# Plotting
plotting(ergmest, vcergmest)

# Hypothesis Test
NBoot = 1000

# Estimate under H0
est0 = estimate_vcergm(object = object, network = networks, attr = attr,
                       degree.spline = degree.spline, interior.knot = interior.knot, 
                       directed = directed, constant = TRUE) # Constant = TRUE (No temporal heterogeneity)

# Calculate test statistic
teststat = test_statistic(object = object, networks = networks, attr = attr, 
                          phi0 = est0$phi.hat, phi1 = vcergmest$phi.hat, 
                          directed = directed)

# P-value
boot = bootstrap_test(object = object, networks = networks, phicoef0 = est0$phicoef.hat, teststat = teststat,
                      degree.spline = degree.spline, interior.knot = interior.knot, directed = directed, 
                      NBoot = NBoot)
