# setwd("~/GRPC")
setwd("C://Users//Liam//Documents//GRPC")

source("Scripts//Init_GRPC.R", echo=TRUE)

jags.data <- list(N.YEARS = N.YEARS,
				  LAMBDA.VR = rep(LAMBDA.VR, times=N.YEARS),
				  COUNTS = COUNTS,
				  P.SURVEYED = P.SURVEYED)

out.params <- c("a",
				"b",
				"sigma.lambda",
				"lambda.lek")

ipm <- jags(data=jags.data, 
			parameters.to.save=out.params,
			model.file="Scripts//ipm.GRPC.jags",
			n.chains=3,
			n.thin=5,
			n.iter=20000,
			n.burnin=10000,
			parallel=TRUE)
