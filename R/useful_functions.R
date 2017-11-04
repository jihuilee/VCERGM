
# Miscellaneous functions

# Calculating quantiles  
quan.data.frame = function(dat, stat, method) # quantiles = c(0.25, 0.5, 0.75)
{
  K = nrow(dat)
  missing.indx = which(is.na(dat[,1]))
  available.indx = setdiff(1:K, missing.indx)
  quan.fun = function(x){quantile(as.numeric(x), probs = c(0.25, 0.5, 0.75), na.rm = TRUE)}
  quan.value = apply(dat, 1, quan.fun)
  out = data.frame(Q1 = quan.value[1,], Median = quan.value[2,], Q3 = quan.value[3,],
                   Time = 1:K, Stat = rep(stat, K), Method = rep(method, K))
  return(out)
}

# Extract legend
g_legend = function(a.gplot) {
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)}
