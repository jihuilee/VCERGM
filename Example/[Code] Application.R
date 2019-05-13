#################################################
## Application (Voting data)                    #
## Calculate the estimates of                   #
##         a) Cross-sectional ERGMs             #
##         b) ad hoc 2-step procedure           #
##         c) TERGM                             #
##         d) VCERGM                            #
#################################################

rm(list=ls())

library(VCERGM)

# Load the dataset
# U.S. Congress Political Voting Data
# Load the data (Congress 40 - 113)
# load("Rollcall.RData")

# Load the required packages
library(splines)
library(ergm)
library(ggplot2)
library(gridExtra)
library(statnet.common)
library(btergm)

# Dataset
networks = Rollcall$networks # Networks
attr = Rollcall$attr # Political affiliation

# Model setup: 
# Network statistics used for fitting models are 1) Triangle, 1) Two-star, and 3) Political affiliation (nodemix)
object = Net ~ triangle + kstar(2) + nodemix("attr1")
stat = c("Triangle", "Two-star", "Dem & Dem", "Dem & Rep", "Rep & Rep")
directed = FALSE # Undirected network

# Degree of spline and number of knots for basis expansion
degree.spline = 3
interior.knot = 30

# Calculate the estimates of       
# a) Cross-sectional ERGMs       
# b) ad hoc 2-step procedure           
# c) VCERGM     
# d) TERGM

#####################################################
# Cross-sectional ERGMs and ad hoc 2-step procedure #
#####################################################

ergmest = cross_sectional_ergm(object = object, network = networks, attr = attr, directed = directed, 
                               degree.spline = degree.spline, interior.knot = interior.knot)

crossERGM_est = ergmest$phi.hat # Cross-sectional ERGMs: hat(phi(t))
adhoc_est = ergmest$phi.hat.smooth # Ad hoc 2-step procedure: hat(phi(t))

##########
# VCERGM #
########## 

vcergmest = estimate_vcergm(object = object, network = networks, attr = attr,
                            degree.spline = degree.spline, interior.knot = interior.knot, 
                            directed = directed, constant = FALSE)

vcergm_est = vcergmest$phi.hat # VCERGM: hat(phi(t))

######### 
# TERGM #
######### 

# Convert a list of adjacency matrices into a list of network objects
networks2 = from_adj_to_net(networks, attr, directed)
object2 = nonsimp_update.formula(object, networks2 ~ ., from.new = TRUE)

tergmest = btergm(object2, R = 200)

tergm_est = attr(tergmest, "coef") # TERGM: hat(phi)

#nstat = nrow(ergm.phi.hat) # Number of network statistics
stat = rownames(ergm.phi.hat) # Name of network statistics
#ntime = length(timeseq) # Number of time points

# Plotting the results
plotting2(ergmest, vcergmest, tergmest, label = stat, timeseq = 40:113, xlab = "Congress")