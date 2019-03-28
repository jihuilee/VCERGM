#' Calculate network statistics for temporal networks
#'
#' @param networks A list of observed networks. It should have a list() object.
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param attr Node-specific attributes. Default is NULL.
#' 
#' @importFrom network network
#' @importFrom statnet.common nonsimp_update.formula
#' @export

calculate_netstat = function(networks, object, attr = NULL)
{
  netstat = NULL
  for (s in 1:length(networks)) {
    nets = networks[[s]]
    if (is.null(attr) == FALSE) {
      if (is.vector(attr[[s]])) {
        attr.s = vector("list", 1)
        attr.s[[1]] = attr[[s]]
        names(attr.s) = "attr1"
      } else {
        attr.s = vector("list", ncol(attr[[s]]))
        for (l in 1:ncol(attr[[s]])) {
          attr.s[[l]] = attr[[s]][, l]
        }
        names(attr.s) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
      nets = network(nets, vertex.attr = attr.s, directed = directed)
    } else {
      nets = network(nets, directed = directed)
    }
    z = deparse(object[[3]])
    formula.s = nonsimp_update.formula(object, nets ~ ., from.new = TRUE)
    temp = summary(formula.s)
    netstat = cbind(netstat, temp)
  }
  colnames(netstat) = 1:length(networks)
  return(netstat)
}

