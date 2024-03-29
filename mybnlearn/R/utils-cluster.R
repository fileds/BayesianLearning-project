# check whether the cluster is running.
isClusterRunning = function(cl) {

  tryCatch(any(unlist(parallel::clusterEvalQ(cl, TRUE))),
    error = function(err) { FALSE })

}#ISCLUSTERRUNNING

# check the status of the snow/parallel cluster.
check.cluster = function(cluster) {

  if (missing(cluster) || is.null(cluster))
    return(NULL)
  if (!is(cluster, "cluster"))
    stop("cluster is not a valid cluster object.")
  check.and.load.package("parallel")
  if (!isClusterRunning(cluster))
    stop("the cluster is stopped.")

  return(cluster)

}#CHECK.CLUSTER

# get the number of slaves.
nSlaves = function(cluster) {

  length(cluster)

}#NSLAVES

slaves.setup = function(cluster) {

  # set the test counter in all the cluster nodes.
  parallel::clusterEvalQ(cluster, library(bnlearn))
  parallel::clusterEvalQ(cluster, reset.test.counter())

}#SLAVE.SETUP

# smart parSapply() that falls back to standard sapply(), but with defaults to
# simplify = FALSE.
smartSapply = function(cl, ..., simplify = FALSE, USE.NAMES = TRUE) {

  if (is.null(cl))
    sapply(..., simplify = simplify, USE.NAMES = USE.NAMES)
  else
    parallel::parSapplyLB(cl = cl, ..., simplify = simplify, USE.NAMES = USE.NAMES)

}#SMARTSAPPLY

