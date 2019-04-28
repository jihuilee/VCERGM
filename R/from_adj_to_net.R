#' Convert a list of adjacency matrices into a list of network objects
#'
#' @param networks A list of adjacency matrices. It should have a list() object.
#' @param attr attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#'
#' @importFrom network network
#' @export

from_adj_to_net = function(networks, attr = NULL, directed = FALSE)
{
  out = vector("list", length(networks))

  for (s in 1:length(networks)) {

    nets = networks[[s]]

    if (is.null(attr) == FALSE)
    {
      if (is.vector(attr[[s]])) {
        attr.s = vector("list", 1); attr.s[[1]] = attr[[s]]; names(attr.s) = "attr1"
      } else {attr.s = vector("list", ncol(attr[[s]]))
      for (l in 1:ncol(attr[[s]])) {attr.s[[l]] = attr[[s]][,l]}
      names(attr.s) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
      nets = network(nets, vertex.attr = attr.s, directed = directed)
    } else{nets = network(nets, directed = directed)}

    out[[s]] = nets
  }
  return(out)
}
