#' Plot the (ERGM, VCERGM, TERGM) results
#'
#' @param ergmest ERGM estimate object
#' @param vcergmest VCERGM estimate object
#' @param tergmest TERGM estimate object
#' @param label Label(s) for plot title. Default is NULL (Use the label from VCERGM object).
#' @param interval Interval size for time in x-axis. Default is 10.
#' @param timeseq User defined timeseq. Default is NULL (Time starts from 1).
#' @param xlab Name for x-axis. Default is "Time".
#' @param theme ggplot theme. Default is NULL.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot_build
#' @importFrom gridExtra grid.arrange
#' @importFrom gridExtra arrangeGrob
#' @export

plotting2 = function(ergmest, vcergmest, tergmest, interval = 10, xlab = "Time", label = NULL, timeseq = NULL, theme = NULL)
{
  ergm.phi.hat = ergmest$phi.hat
  ergm.phi.hat.smooth = ergmest$phi.hat.smooth
  vcergm.phi.hat = vcergmest$phi.hat
  tergm.phi.hat = attr(tergmest, "coef")

  # Number of time points
  if(is.null(timeseq)){timeseq = 1:ncol(ergm.phi.hat)}

  nstat = nrow(ergm.phi.hat) # Number of network statistics
  stat = rownames(ergm.phi.hat) # Name of network statistics
  ntime = length(timeseq) # Number of time points

#  if (nstat > 3) {NRow = 2} else {NRow = 1}

  plot.dat = data.frame(Value = c(c(t(ergm.phi.hat)), c(t(ergm.phi.hat.smooth)), c(t(vcergm.phi.hat)), rep(tergm.phi.hat, each = ntime)),
                        Method = rep(c("ERGM", "ERGM2", "VCERGM", "TERGM"),
                                       each = nstat * ntime),
                        Time = rep(timeseq, nstat * 4),
                        Stat = rep(rep(stat, each = ntime), 4))

  plot.dat$Method = factor(plot.dat$Method, levels = c("ERGM", "ERGM2", "VCERGM", "TERGM"))
  plist = vector("list", length(stat))

  for (i in 1:nstat)
  {
    plot.dat.i = plot.dat[plot.dat$Stat == stat[i], ]
    plist[[i]] = ggplot(data = plot.dat.i, aes(x = as.factor(Time), y = Value, col = Method, group = Method)) +
                  geom_line() + geom_point(data = plot.dat.i[plot.dat.i$Method != "TERGM", ], size = 1) + xlab("Time") + ylab(expression(hat(phi)(t))) +
                  scale_x_discrete(breaks = c(1, interval * (1:floor(max(timeseq)/interval)))) +
                  ggtitle(stat[i]) + theme(legend.position = "bottom") + xlab(xlab)

    if(!is.null(label)) {plist[[i]] = plist[[i]] + ggtitle(label[i])}
    if(!is.null(theme)) {plist[[i]] = plist[[i]] + theme}
  }

  g_legend = function(a.gplot) {
    tmp = ggplot_gtable(ggplot_build(a.gplot))
    leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend = tmp$grobs[[leg]]
    return(legend)}

  mylegend = g_legend(plist[[1]])

  for (i in 1:nstat) {plist[[i]] = plist[[i]] + theme(legend.position = "none") +
                                    theme(plot.title = element_text(hjust = 0.5))}

  plots = do.call("arrangeGrob", c(plist, ncol = nstat))
  print(grid.arrange(plots, mylegend, heights = c(9/10, 1/10)))

#  plot2 = ggplot(data = plot.dat, aes(x = as.factor(Time), y = Value,
#                      col = Method, group = Method), alpha = 0.8) +
#                      geom_line() + geom_point(size = 1) +
#                      facet_grid(. ~ Stat) + xlab("Time") + ylab("Phi")
#  print(plot2)
}
