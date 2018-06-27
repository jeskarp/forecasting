#' Connect to the dide hpc cluster
#'
#' This function accepts the necessary arguments to connect to the dide cluster
#' and returns a queue object, which you can use to submit tasks, query the
#' status and collect results.
#'
#' @param dir the directory (on your Q: drive) in which the queue is set up
#'
#' @param sources the path to the file containing your functions
#'
#' @param pkgs a list of the package names you use
#'
#' @param pkg_sources an optional \code{package_sources} object from the package
#'   \code{provisionr}, describing where to install non-CRAN packages from
#'
#' @param contexts the name of the directory the queue and outputs will be saved in
#'
#' @param clust the shorthand name for the cluster you want to connect to - must
#'   match to a name in the \code{clusters} object defined within the function
#'
#' @seealso \link{https://mrc-ide.github.io/didehpc/vignettes/didehpc.html} for
#'   full documentation of the packages
#' 
#' @examples
#'
#' dir <- "Q://my_cluster_directory"
#' sources <- "my_functions.R"
#' pkgs <- c("ggplot2", "outbreaker2", "incidence")
#' pkg_sources <- provisionr::package_sources(local=c("~/outbreaker2_1.1.0.tar.gz"),
#' github = c("reconhub/incidence")
#'
#' obj <- cluster_connect(dir, sources, pkgs, pkg_sources, 'queue_1')
#'
#' ## Send off a single run of the function "my_func"
#' t <- obj$enqueue(my_func(arg1, arg2))
#'
#' ## Send off multiple runs
#' par <- data.frame(arg1 = runif(10, 0, 1), arg2 = runif(10, 1, 2))
#' q <- obj$enqueue_bulk(par, my_func)
#' 
cluster_connect <- function(dir = "/didetmp/Janetta/", sources, pkgs, pkg_sources = NULL, contexts = 'test',
                            clust = 'hn') {

  if(!require(didehpc)) source("https://mrc-ide.github.io/didehpc/install")
  
  setwd(dir)

  clusters <- c(hn = "fi--dideclusthn",
                mrc = "fi--didemrchnb")
  
  options(didehpc.home = dir,
          didehpc.cluster = clusters[clust],
          didehpc.credentials = "~/.smbcredentials")

  didehpc::web_login()

  ctx <- context::context_save(contexts, packages = pkgs,
                               sources = sources,
                               package_sources = pkg_sources)

  obj <- didehpc::queue_didehpc(ctx)

  return(obj)

}
