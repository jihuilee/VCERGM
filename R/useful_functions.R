
# Miscellaneous functions

# Calculating quantiles  
quan.data.frame = function(dat, stat, method) # quantiles = c(0.25, 0.5, 0.75)
{
  T = nrow(dat)
  missing.indx = which(is.na(dat[,1]))
  available.indx = setdiff(1:T, missing.indx)
  quan.fun = function(x){quantile(as.numeric(x), probs = c(0.25, 0.5, 0.75), na.rm = TRUE)}
  quan.value = apply(dat, 1, quan.fun)
  out = data.frame(Q1 = quan.value[1,], Median = quan.value[2,], Q3 = quan.value[3,],
                   Time = 1:T, Stat = rep(stat, T), Method = rep(method, T))
  return(out)
}

# Extract legend
g_legend = function(a.gplot) {
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)}

# (Transform) Correlation & Precision
transform = function(dat, window = 200, overlap = 100)
{
  timepoint = (nrow(dat) - overlap)/(window - overlap)
  
  Correlation = Precision = vector("list", timepoint)
  for (i in 1:timepoint)
  {
    loc1 = (i - 1) * (window - overlap) + 1
    loc2 = i * window - (i - 1) * overlap
    
    dat.i = dat[loc1:loc2, ]
    
    corr.i = cor(dat.i)
    diag(corr.i) = 0
    Correlation[[i]] = corr.i
    
    prec.i = solve(cov(dat.i)) # standardize it to make it partial
    diag(prec.i) = 0
    Precision[[i]] = prec.i
  }
  return(list(Correlation = Correlation, Precision = Precision))
}

# (Binary networks) Dichotomize correlation / precision
threshold = function(networks, density = 0.1)
{
  new = vector("list", length(networks))
  for(i in 1:length(networks))
  {
    net.i = networks[[i]]
    thres = as.numeric(quantile(net.i, probs = 1 - density))
    new[[i]] = (net.i > thres) + 0
  }
  return(new)
}
