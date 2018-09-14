setwd("~/GRPC")

require(tidyverse)
require(lubridate)
require(magrittr)

# --------------------------------
# Initialize matrix vital rates
# from "random LTRE LL gpc.r"
# SM poplation only
# --------------------------------

# the original list, from "random LTRE LL gpc.r" :
# SM_vitalrates<- list(
#   SMnest = 0.86
#   SMrenest  =  0.4,
#   SMtcl1    =  13.4,
#   SMtcl2    =  10.8,
#   SMnsurv1y =  0.16,
#   SMnsurv2y =  0.28,
#   SMnsurv1a =  0.18,
#   SMnsurv2a =  0.31,
#   SMcpere   =  0.83,
#   SMbsurv  =   0.34,
#   SMfperc   =  0.51,
#   SMsj      =  0.38,
#   SMsy      =  0.34,
#   SMsa      =  0.42)

# NOTE N.Prob (old: SMnest) and Re.N.Prob (old: SMrenest) are here set to 1

N.Prob = 1		# nesting probability per female
Re.N.Prob = 1	# re-nesting probability per female

# N.Prob = 0.86		# nesting probability per female
# Re.N.Prob = 0.4		# re-nesting probability per female

C.1.Size = 13.4		# clutch size (#eggs) per first nest
C.2.Size = 10.8		# clutch size (#eggs) per second nest

N.Surv.1.Y = 0.16 	# nest survival probability for yearlings' first nest
N.Surv.2.Y = 0.28	# nest survival probability for yearlings' second nest

N.Surv.1.A = 0.18	# nest survival probability for adults' first nest
N.Surv.2.A = 0.31	# nest survival probability for adults' second nest

B.Surv = 0.34		# brood survival probability

HS = 0.83			# hatching success (Chicks / Eggs)
FS = 0.51			# fledging success (Fledglings / Chicks)

S.J = 0.38			# Juvenile Female survival
S.Y = 0.34			# Yearling Female survival
S.A = 0.42			# Adult Female survival

# Fecundity of yearlings
F.Y = ((N.Prob * C.1.Size * N.Surv.1.Y) +
		(1 - N.Surv.1.Y) *
		(Re.N.Prob * C.2.Size * N.Surv.2.Y)) *
	  HS * B.Surv * FS * 0.5

# Fecundity of adults
F.A = ((N.Prob * C.1.Size * N.Surv.1.A) +
		 (1 - N.Surv.1.A) *
		 (Re.N.Prob * C.2.Size * N.Surv.2.A)) *
	  HS * B.Surv * FS * 0.5

# Translate to matrix to determine eigenvalues->lambda
values.MAT <- list(F.Y = F.Y,
					F.A = F.A,
					S.J = S.J,
					S.Y = S.Y,
					S.A = S.A)

elements.MAT <- expression(F.Y * S.J, F.A * S.J,
						    S.Y   ,    S.A)

Vrs.MAT <- sapply(elements.MAT, eval, values.MAT, NULL)

A.MAT <- matrix(Vrs.MAT, nrow=2, byrow=TRUE)

evs <- eigen(A.MAT)

# LAMBDA
imax <- which.max(evs$values)
lambda <- Re(evs$values[imax])


# --------------------------------
# Initialize m
# from "UK_lek_metadata_tables.xlsx" original datasheet
# --------------------------------

leks <- read.csv("Data/lek_metadata.csv",
				 col.names=c("trapped",
				 			 "edited",
				 			 "id",
				 			 "date",
				 			 "temp",
				 			 "wind",
				 			 "sky",
				 			 "n.tot",
				 			 "n.m",
				 			 "n.f"),
				 colClasses=c("logical", 
				 			  "character",
				 			  "character", 
				 			  "character",
				 			  "integer",
				 			  "integer",
				 			  "factor",
				 			  "character",
				 			  "character",
				 			  "character"),
				 stringsAsFactors=FALSE)

leks$date <- mdy(leks$date)
leks$n.tot[leks$n.tot=="#VALUE!"] <- NA
leks$n.m[leks$n.m=="#VALUE!"] <- NA
leks$n.f[leks$n.f=="#VALUE!"] <- NA

leks$n.tot <- as.numeric(leks$n.tot)
leks$n.m <- as.numeric(leks$n.m)
leks$n.f <- as.numeric(leks$n.f)
# Lek count data (follow Winder et al. 2015)
# Applied to both sexes

# (a) find the mean trap vs. flush coefficient for total pop sizes
trap_v_flush <- leks %>%
					group_by(trapped, id, year(date)) %>%
					summarise(max.m = max(n.m, na.rm=TRUE),
							  max.f = max(n.f, na.rm=TRUE))

males <- leks %>%
			group_by(id, year(date)) %>%
			summarise(max.m = max(n.m, na.rm=TRUE))

counts <- leks %>%
			mutate(n.tot = as.numeric(n.tot),
				   n.m = as.numeric(n.m),
				   n.f = as.numeric(n.f),
				   date = mdy(date)) %>%
			group_by(id, year(date)) %>%
			summarize(mean.tot=mean(n.tot, na.rm=TRUE),
					  mean.m = mean(n.m, na.rm=TRUE),
					  mean.f = mean(n.f, na.rm=TRUE))


# TODO 
# From winder et al. 2015
# FOR each sex
# Find max number of males during trap vs. flush for each lek per year
# Discount Flush counts to equalize
# calculate weighted mean given number of visits of each kind

