# setwd("~/GRPC")
setwd("C://Users//Liam//Documents//GRPC")

source("Scripts/Init_GRPC.R", echo=TRUE)

numChains <- 3
numThin <- 20
numBurn <- 50000
numIters <- 100000

jags.data <- list(N.YEARS    = N.YEARS,
		  COUNTS     = COUNTS,
		  P.SURVEYED = P.SURVEYED,
		  N.PROB     = N.PROB, 
		  RE.N.PROB  = RE.N.PROB,
		  C.1.SIZE   = C.1.SIZE,
		  C.2.SIZE   = C.2.SIZE,
		  N.SURV.1   = N.SURV.1,
		  N.SURV.2   = N.SURV.2,
		  HS.1       = HS.1,
		  HS.2       = HS.2,
		  B.SURV     = B.SURV,
		  FS         = FS,
		  S.J        = S.J,
		  S.A        = S.A)

jags.inits <- function() {
	list(m.numImms=runif(1,0,100),
		 N.Imm=c(NA,round(runif(N.YEARS-1,0,100))))
}

out.params <- c("f",
		"N.a",
		"N.i",
		"N.tot",
		"omega",
		"lambda")

ipm <- jags(data=jags.data,
	    inits=jags.inits,
	    parameters.to.save=out.params,
	    model.file="Scripts//ipm_grpc.jags",
	    n.chains=numChains,
	    n.thin=numThin,
	    n.burnin=numBurn,
	    n.iter=numIters,
	    parallel=TRUE)

# summarize ipm output for analyses/visualization
out <- as.data.frame(ipm$summary)
out <- bind_cols(as.data.frame(rownames(out)), out)
colnames(out)[1] <- "param"

out <- separate(out, param, c("param", "year"), sep="\\[|\\]", 
                extra="drop", fill="right")

colnames(out) <- c("param", "year", "mean", "sd", "low", 
                   "twentyfive", "mid", "seventyfive", "high", 
                   "rhat", "neff", "overlap", "f")

out$year <- as.numeric(out$year) + 2006

ipm_sim <- ipm$sims.list    # the raw results at every iterations
