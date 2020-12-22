#!/usr/bin/env Rscript
library(foreach)
library(parallel)
library(doParallel)
library(conStruct)


kdm.x.validation <- function(train.prop = 0.9, n.reps, K, freqs = NULL, data.partitions = NULL, geoDist, coords, prefix, n.iter, make.figs = FALSE, save.files = FALSE, n.nodes=16, parallel=T, ...) {
    if (is.null(n.nodes))
        n.nodes = as.integer(Sys.getenv("PBS_NCPUS", parallel::detectCores(logical=F)))
    call.check <- conStruct:::check.xval.call(args <- as.list(environment()))
    data.partitions <- conStruct:::make.data.partitions(n.reps,freqs,train.prop)
    conStruct:::check.data.partitions.arg(args <- as.list(environment()))
    cat("Data partitions created successfully\n")
    #.par.clust = makeCluster(n.nodes,type="FORK")
    #registerDoParallel(.par.clust)
    n.mdl = 2 * n.reps * length(K)
    cat(paste("Using", n.nodes, "cores to run", n.mdl, "models\n"))
    x.val = foreach::foreach(rep.no=1:n.reps) %:%
        foreach::foreach(k=K) %:%
        foreach::foreach(mdl=c("sp", "nsp")) %dopar% {
        spatial = mdl == "sp"
        cs = conStruct:::xval.conStruct(spatial = TRUE, K = k, 
                        data = data.partitions[[rep.no]]$training, 
                        geoDist = geoDist, coords = coords, 
                        prefix = paste0(prefix, "_", mdl, "_", "rep", rep.no, "K", k), 
                        n.iter = n.iter, make.figs = make.figs, save.files = save.files)
        cat(paste0("mdl=", mdl, " K=",k, " rep=", rep.no, " Done\n"))
        list(mdl=mdl, K=k, rep=rep.no, construct=cs)
    }
    save(data.partitions,file=paste0(prefix, ".xval.data.partitions.Robj"))
    save(x.val,file=paste0(prefix, ".xval.results.Robj"))
    #stopCluster(.par.clust)
    x.val
}

