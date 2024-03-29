
# operate on arcs in a network.
arc.operations = function(x, from, to, op = NULL, check.cycles,
    check.illegal, update = TRUE, debug = FALSE) {

  available.ops = c("set", "drop", "reverse", "seted", "droped")

  # check x's class.
  check.bn(x)
  # check the op code.
  if (op %!in% available.ops)
    stop("valid op codes are 'set', 'drop' and 'reverse'.")
  # a valid node is needed.
  check.nodes(nodes = from, graph = x, max.nodes = 1)
  # another valid node is needed.
  check.nodes(nodes = to, graph = x, max.nodes = 1)
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")
  # check logical flags (debug, check.cycles, update).
  check.logical(debug)
  check.logical(check.cycles)
  check.logical(update)

  # add/reverse/orient the arc (or the edge).
  if (op == "set") {

    if (debug)
      cat("* setting arc", from, "->", to, ".\n")
    if (check.illegal && is.listed(x$learning$illegal, c(from, to)))
      stop("arc ", from, " -> ", to,
        " is not valid due to the parametric assumptions of the network.")

    x$arcs = set.arc.direction(from, to, x$arcs, debug = debug)

  }#THEN
  else if (op == "drop") {

    if (debug)
      cat("* dropping any arc between ", from, "and", to, ".\n")

    x$arcs = drop.arc.backend(x$arcs, c(from, to), debug = debug)

  }#THEN
  else if (op == "reverse") {

    if (debug)
      cat("* reversing any arc between ", from, "and", to, ".\n")
    if (check.illegal && is.listed(x$learning$illegal, c(to, from)))
      stop("arc ", to, " -> ", from,
        " is not valid due to the parametric assumptions of the network.")

    x$arcs = reverse.arc.backend(from, to, x$arcs, debug = debug)

  }#THEN
  else if (op == "seted") {

    if (debug)
      cat("* setting undirected arc", from, "-", to, ".\n")
    if (check.illegal && is.listed(x$learning$illegal, c(to, from), either = TRUE))
      stop("undirected arc ", to, " - ", from,
        " is not valid due to the parametric assumptions of the network.")

    x$arcs = set.edge.backend(from, to, x$arcs, debug = debug)

  }#THEN
  else if (op == "droped") {

    if (debug)
      cat("* dropping undirected arc", from, "-", to, ".\n")

    x$arcs = drop.edge.backend(x$arcs, c(from, to), debug = debug)

  }#THEN

  # check whether the graph contains directed cycles; not needed if an arc
  # is dropped.
  if (check.cycles && (op != "drop"))
    if (!is.acyclic(x$arcs, names(x$nodes), debug = debug, directed = TRUE))
      stop("the resulting graph contains cycles.")

  # update the network structure.
  if (update) {

    # build the adjacency matrix only once.
    amat = arcs2amat(x$arcs, names(x$nodes))
    # check which nodes have to be updated.
    updated.nodes = unique(c(from, to, x$nodes[[from]]$mb, x$nodes[[to]]$mb))
    # update the chosen nodes.
    for (node in updated.nodes)
      x$nodes[[node]] = cache.partial.structure(names(x$nodes),
        target = node, amat = amat, debug = debug)

  }#THEN

  invisible(x)

}#ARC.OPERATIONS

