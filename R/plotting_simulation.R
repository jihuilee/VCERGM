#' Plot the ERGM and VCERGM results when nsim > 1
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'. 
#' @param true.phi TRUE phi(t) functions (number of network statistics x K vector) 
#' @param ergm.phi.hat Estimated phi(t) from cross-sectional ERGMs. It should be a list of length nsim. For each simulation, the estimated phi.hat is a (number of network statistics x K) matrix.
#' @param ergm.phi.hat.smooth Estimated phi(t) from ad hoc 2-step procedure. 
#' @param vcergm.phi.hat Estimated phi(t) from VCERGM.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom gridExtra arrangeGrob
#' @export

plotting_simulation = function(object, true.phi, ergm.phi.hat, ergm.phi.hat.smooth, vcergm.phi.hat)
{
# Number of simulations & time points
nsim = length(ergm.phi.hat)
K = ncol(ergm.phi.hat[[1]])

# Network statistics
stat = unlist(strsplit(deparse(object[[3]]), " "))
stat = stat[!stat %in% c("+", "=", "TRUE)", "FALSE)")]
nstat = length(stat)

Summary = matrix(NA, nrow = nstat, ncol = 3) # Compare mean/sd of ERGM, ERGM2, and VCERGM

# True
plot.true = data.frame(Value = c(t(true.phi)),
                       Time = rep(1:K, nstat),
                       Stat = rep(stat, each = K),
                       Method = rep("True", nstat * K))

# Estiamtes  
plot.frame = data.frame(Stat = rep(stat, each = K), Time = rep(1:K, nstat))
plot.ergm0 = plot.ergm.smooth0 = plot.vcergm0 = rep(list(matrix(NA, nrow = K, ncol = nsim)),nstat)
plot.ergm = plot.ergm.smooth = plot.vcergm = plot.frame

# Concatenate all estimates
for (i in 1:nsim) 
{
  ergm.i = c(t(ergm.phi.hat[[i]]))
  ergm.smooth.i = c(t(ergm.phi.hat.smooth[[i]]))
  vcergm.i = c(t(vcergm.phi.hat[[i]]))
  
  plot.ergm = cbind(plot.ergm, ergm.i)
  plot.ergm.smooth = cbind(plot.ergm.smooth, ergm.smooth.i)
  plot.vcergm = cbind(plot.vcergm, vcergm.i)
  
  for (p in 1:nstat)
  {
    plot.ergm0[[p]][,i] = result$ergm.phi.hat[[i]][p,]
    plot.ergm.smooth0[[p]][,i]  = result$ergm.phi.hat.smooth[[i]][p,]
    plot.vcergm0[[p]][,i]  = result$vcergm.phi.hat[[i]][p,]
  }
  
  index = rep(i, K * nstat)
  ergm.i2 = cbind(plot.frame, ergm.i, index)
  vcergm.i2 = cbind(plot.frame, vcergm.i, index)
}

pplist1 = pplist2 = vector("list", nstat)

for (u in 1:nstat)
{
  stat.u = stat[u]
  ergm.temp0 = plot.ergm0[[u]]
  ergm.temp = plot.ergm.smooth0[[u]]
  vcergm.temp = plot.vcergm0[[u]]
  
  true.u = get(paste("true.phi", u, sep = ""))
  
  # IAE
  ergm.bias0 = apply(ergm.temp0, 2, function(x){sum(abs(x - true.u), na.rm = T)})
  ergm.bias = apply(ergm.temp, 2, function(x){sum(abs(x - true.u))})
  vcergm.bias = apply(vcergm.temp, 2, function(x){sum(abs(x - true.u))})
  
  # Mean / sd
  Summary[u, 1] = paste(round(mean(ergm.bias0), 2), " (", round(sd(ergm.bias0), 2), ")", sep = "")
  Summary[u, 2] = paste(round(mean(ergm.bias), 2), " (", round(sd(ergm.bias), 2), ")", sep = "")
  Summary[u, 3] = paste(round(mean(vcergm.bias), 2), " (", round(sd(vcergm.bias), 2), ")", sep = "")
}

plist = vector("list", nstat)

for (k in 1:nstat)
{
  stat.k = stat[k]
  dat.ergm = plot.ergm0[[k]]
  dat.ergm.smooth = plot.ergm.smooth0[[k]]
  dat.vcergm = plot.vcergm[plot.vcergm$Stat == stat.k, -c(1:2)]
  
  true.k = plot.true[plot.true$Stat == stat.k, ]
  ergm.k = quan.data.frame(dat.ergm, stat.k, "ERGM")
  ergm.smooth.k = quan.data.frame(dat.ergm.smooth, stat.k, "ERGM2")
  vcergm.k = quan.data.frame(dat.vcergm, stat.k, "VCERGM")
  
  quan.dat = rbind(ergm.k, ergm.smooth.k, vcergm.k)
  plist[[k]] = ggplot() + geom_ribbon(data = quan.dat, aes(x = as.factor(Time), ymin = Q1, ymax = Q3, 
                                                           group = Method, col = Method, fill = Method), alpha = 0.3, col = NA) +
    scale_x_discrete(breaks = c(1, 10 * (1:floor(T/10)), T)) +
    geom_line(data = quan.dat, aes(x = as.factor(Time), y = Median, group = Method, col = Method)) +
    geom_line(data = true.k, size = 0.5, alpha = 0.7, aes(x = as.factor(Time), y = Value, group = Method)) +
    xlab("Time") + ylab("") + theme(legend.position = "bottom")
  plist[[k]] = plist[[k]] + ggtitle(stat[k]) + theme(plot.title = element_text(hjust = 0.5))
}

colnames(Summary) = c("ERGM", "ERGM2", "VCERGM")
rownames(Summary) = stat

mylegend = g_legend(plist[[1]])
for (i in 1:length(plist)) {plist[[i]] = plist[[i]] + theme(legend.position = "none")}

plots = do.call("arrangeGrob", c(plist, ncol = nstat))
grid.arrange(plots, mylegend, heights = c(9/10, 1/10))

return(Summary = Summary)
}
