#' Goodness-of-fit diagnostics on VCERGM
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param object2 A formula object of the form (network object) ~ <model terms> for GoF. By default, it is the same as object.
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#' @param netstat A vector of GoF network statistics names. Default is the name from ERGM estimates.
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param interior.knot Number of interior knots to create splines. Default is 10.
#' @param log Logarithm for observed and simulated network statistics. Default is FALSE.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 facet_wrap
#' @export

gof_vcergm = function(object, object2 = NULL, networks, attr = NULL, directed = FALSE, netstat = NULL,
                      degree.spline = 3, interior.knot = 10, nsim = 100, seed = 1234, log = FALSE)
{

  if(is.null(object2)){object2 = object}

  # VCERGM estimates
  vcergmest = estimate_vcergm(object = object, network = networks, attr = attr,
                              degree.spline = degree.spline, interior.knot = interior.knot,
                              directed = directed, constant = FALSE)

  # Simulate networks
  vcergm_phi = vcergmest$phi.hat
  sim.net = simulate_vcergm(object = object, attr = attr, num.nodes.K = unlist(lapply(networks, nrow)), phi = vcergm_phi,
                            nsim = nsim, directed = directed, seed = seed)$Networks

  # Network statistics: observed vs simulated
  obs.netstat0 = calculate_netstat(networks, object2, attr = attr)
  if(is.null(netstat)){netstat = rownames(obs.netstat0)}

  obs.netstat = data.frame(value = c(t(obs.netstat0)), time = rep(1:length(networks), length(netstat)),
                           stat = rep(netstat, each = length(networks)))
  obs.netstat$stat = factor(obs.netstat$stat, labels = netstat)

  sim.netstat0 = lapply(sim.net, function(x){calculate_netstat(networks = x, object = object2, attr = attr)})
  sim.netstat = NULL
  for(t in 1:length(sim.netstat0))
  {
    dat.t = data.frame(value = c(t(sim.netstat0[[t]])), simnum = t, time = rep(1:length(networks), length(netstat)),
                       stat = rep(netstat, each = length(networks)))
    sim.netstat = rbind(sim.netstat, dat.t)
  }

  if(log){
    obs.netstat$value = log(obs.netstat$value)
    sim.netstat$value = log(sim.netstat$value)
  }

  out = ggplot() + geom_boxplot(data = sim.netstat, aes(x = as.factor(time), y = value, fill = stat), alpha = 0.1) +
    theme_bw() + theme(legend.position = "none") + labs(x = "Time", y = "Network statistics") +
    scale_x_discrete(breaks = c(1, 10 * (1:floor(length(networks)/10)), length(networks))) +
    geom_line(data = obs.netstat, aes(x = time, y = value, color = stat)) + facet_wrap(.~stat, scales = "free_y")
  return(out)
}
