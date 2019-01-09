# setwd("~/GRPC")
setwd("C://Users//Liam//Documents//GRPC")

source("Scripts/Init_GRPC.R", echo=TRUE)

numChains <- 3
numThin <- 20
numBurn <- 50000
numIters <- 100000

jags.data <- list(N.YEARS = N.YEARS,
				  COUNTS = COUNTS,
				  P.SURVEYED = P.SURVEYED,
				  N.PROB = N.PROB, RE.N.PROB = RE.N.PROB,
				  mu.C.1.SIZE.Y = mu.C.1.SIZE.Y, tau.C.1.SIZE.Y = tau.C.1.SIZE.Y,
				  mu.C.1.SIZE.A = mu.C.1.SIZE.A, tau.C.1.SIZE.A = tau.C.1.SIZE.A,
				  mu.C.2.SIZE.Y = mu.C.2.SIZE.Y, tau.C.2.SIZE.Y = tau.C.2.SIZE.Y,
				  mu.C.2.SIZE.A = mu.C.2.SIZE.A, tau.C.2.SIZE.A = tau.C.2.SIZE.A,
				  mu.N.SURV.1.Y = mu.N.SURV.1.Y, tau.N.SURV.1.Y = tau.N.SURV.1.Y,
				  mu.N.SURV.1.A = mu.N.SURV.1.A, tau.N.SURV.1.A = tau.N.SURV.1.A,
				  mu.N.SURV.2.Y = mu.N.SURV.2.Y, tau.N.SURV.2.Y = tau.N.SURV.2.Y,
				  mu.N.SURV.2.A = mu.N.SURV.2.A, tau.N.SURV.2.A = tau.N.SURV.2.A,
				  mu.HS = mu.HS, tau.HS = tau.HS,
				  mu.B.SURV.Y = mu.B.SURV.Y, tau.B.SURV.Y = tau.B.SURV.Y,
				  mu.B.SURV.A = mu.B.SURV.A, tau.B.SURV.A = tau.B.SURV.A,
				  mu.FS = mu.FS, tau.FS = tau.FS,
				  mu.S.J = mu.S.J, tau.S.J = tau.S.J,
				  mu.S.Y = mu.S.Y, tau.S.Y = tau.S.Y,
				  mu.S.A = mu.S.A, tau.S.A = tau.S.A)

jags.inits <- function() {
	list(m.numImms=runif(1,0,100),
		 N.Imm=c(NA,round(runif(N.YEARS-1,0,100))))
}

out.params <- c("f.y",
				"f.a",
				"f.i",
				"N.y",
				"N.a",
				"N.i",
				"N.tot",
				"omega",
				"lambda",
				"c.1.size.y",
				"n.surv.1.y",
				"c.2.size.y",
				"n.surv.2.y",
				"hs",
				"b.surv.y",
				"fs")

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
