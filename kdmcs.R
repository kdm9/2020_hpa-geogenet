kdm.x.validation <- function(n.reps, K, train.prop = 0.9, freqs = NULL,
                             data.partitions = NULL, geoDist, coords, prefix,
                             n.iter, make.figs = FALSE, save.files = FALSE,
                             n.nodes=16, parallel=T, ...) {
	#######################################
	#  Check args (from check.xval.call)  #
	#######################################
	args = as.list(environment())
	conStruct:::check.for.files(args)
	conStruct:::check.genetic.data.arg(args)
	args$spatial <- TRUE
	conStruct:::check.geoDist.arg(args)
	conStruct:::check.coords.arg(args)

	#####################
	#  Make partitions  #
	#####################
    if (is.null(data.partitions)) {
        data.partitions = conStruct:::make.data.partitions(n.reps, freqs, train.prop)
    }
	conStruct:::check.data.partitions.arg(args = as.list(environment()))

	###############
	#  Main loop  #
	###############
	x.val = foreach::foreach(rep.no=1:n.reps, .combine=bind_rows) %:%
		foreach::foreach(k=K, .combine=bind_rows) %:%
		foreach::foreach(mdl=c("sp", "nsp"), .combine=bind_rows) %dopar% {
            
            # Prep data structures
            spatial = mdl == "sp"
            dat = data.partitions[[rep.no]]
    
            # Run construct mcmc chain
            cs = conStruct:::xval.conStruct(
                spatial = spatial,
                K = k, 
                data = dat$training, 
                geoDist = geoDist,
                coords = coords, 
                prefix = paste0(prefix, "_", mdl, "_", "rep", rep.no, "K", k), 
                n.iter = n.iter,
                make.figs = make.figs,
                save.files = save.files)

            # calculate fit
            fit = unlist(conStruct:::fit.to.test(dat$testing, cs[[1]]))
            mean.fit = mean(fit, na.rm=T)

            # re-create data block (downstream stuff needs it)
            dblock = conStruct:::xval.make.data.block(k, dat$training, coords, spatial, geoDist)
            dblock = conStruct:::unstandardize.distances(dblock)

            cat(paste0("mdl=", mdl, " K=",k, " rep=", rep.no, " Done\n"))
            tibble_row(mdl=mdl, K=k, rep=rep.no,
                       construct.res=cs, mdl.fit=list(fit), mean.fit=mean.fit,
                       data.part=list(dat), data.block = dblock)
	}
	return(x.val)
}
