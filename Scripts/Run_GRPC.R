# setwd("~/GRPC")
# setwd("C://Users//Liam//Documents//GRPC")

source("Scripts//Init_GRPC.R", echo=TRUE)

model <- 1		# 1 or 2
numChains <- 3
numThin <- 20
numBurn <- 500000
numIters <- 1000000

Model1 <- function() {
	jags.data <- list(N.YEARS = N.YEARS,
					  COUNTS = COUNTS,
					  P.SURVEYED = P.SURVEYED,
					  F.SY = F.Y,
					  F.ASY = F.A,
					  F.Imm = mean(F.Y, F.A),
					  S.J = S.J,
					  S.SY = S.Y,
					  S.ASY = S.A,
					  S.Imm = mean(S.Y, S.A))

	jags.inits <- function() {
		list(m.numImms=runif(1,0,100),
			 N.Imm=c(NA,round(runif(N.YEARS-1,0,100))))
	}

	out.params <- c("N.SY",
					"N.ASY",
					"N.Imm",
					"N.tot",
					"omega",
					"lambda")

	ipm <- jags(data=jags.data, 
				inits=jags.inits,
				parameters.to.save=out.params,
				model.file="Scripts//ipm_1.GRPC.jags",
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

	return(list(out=out, ipm_sim=ipm_sim))
}

Model2 <- function() {
	jags.data <- list(N.YEARS = N.YEARS,
					  LOG.LAMBDA.VR = rep(log(LAMBDA.VR), times=N.YEARS),
					  COUNTS = COUNTS,
					  P.SURVEYED = P.SURVEYED)

	out.params <- c("a",
					"b",
					"sigma.lambda",
					"lambda.lek",
					"theta.vr")

	ipm <- jags(data=jags.data, 
				parameters.to.save=out.params,
				model.file="Scripts//ipm_2.GRPC.jags",
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

	return(list(out=out, ipm_sim=ipm_sim))
}

if (model==1) {
	ipm <- Model1()
	out <- ipm$out
	save(out, file="Output/out_Model_1")
} else if (model==2) {
	ipm <- Model2()
	out <- ipm$out
	save(out, file="Output/out_Model_2")
} else {
	print("Uh oh! Specify model 1 or 2")
}
