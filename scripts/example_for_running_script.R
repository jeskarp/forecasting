source("~/Development/forecasting/cluster_connect.R")
obj <- cluster_connect(dir = "/jes213", sources = "test_script.R", pkgs = NULL, clust = 'hn', contexts = 'try2')

obj$cluster_load()

par <- data.frame(x = 0:10)

t <- obj$enqueue(func(5))
t$status()

q <- obj$enqueue_bulk(par, func)
bundle.name <- 'colourific_mayfly'
bundle <- obj$task_bundle_get(bundle.name)
bundle$status()

t$status()

res <- t$result()
hist(res$output)
