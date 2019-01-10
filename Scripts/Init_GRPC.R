# ****************AGE-CLASS BRANCH**************** 
# where lower-level vital rates are extracted from McNew et al. 2012,
# 	which includes differentiation between yearling and adult vital rates
# 	but only covers 3y pre-construction period (2007-2009)

# setwd("~/GRPC")
# setwd("C://Users//Liam//Documents//GRPC") 

# install.packages(c("tidyverse", "lubridate", "magrittr", "jagsUI"))

require(tidyverse)
require(lubridate)
require(magrittr)
require(jagsUI)

# --------------------------------
# Initialize matrix vital rates
# Smoky Hills population only
# Mean and SE values from methods and Table 1 in McNew et al. 2012
# --------------------------------

# Conversion from 
# Standard Error (as in McNew et al. 2012 Table 1)
# 	to Tau [= 1 / sigma^2]
# for proper parameterization of dnorm priors in JAGS
tau <- function(se, n) 
{ 
	sigma = se * sqrt(n)	# standard deviation
	tau = sigma^-2		# 1 / variance
	return (tau)
}

# NOTE set to 1 following McNew et al. 2012 methods
N.PROB        = 1		# nesting probability per female
RE.N.PROB     = 1	    # re-nesting probability per female

# clutch size (#eggs) per first nest (YEARLINGS)
mu.C.1.SIZE.Y = 13.1
se.C.1.SIZE.Y = 0.33
n.C.1.SIZE.Y  = 33

tau.C.1.SIZE.Y = tau(se.C.1.SIZE.Y, n.C.1.SIZE.Y) 

# clutch size (#eggs) per first nest (ADULTS)
mu.C.1.SIZE.A = 13.7
se.C.1.SIZE.A = 0.33
n.C.1.SIZE.A  = 32

tau.C.1.SIZE.A = tau(se.C.1.SIZE.A, n.C.1.SIZE.A)

# clutch size (#eggs) per second nest (YEARLINGS)
mu.C.2.SIZE.Y = 11.2		
se.C.2.SIZE.Y = 0.56
n.C.2.SIZE.Y  = 11

tau.C.2.SIZE.Y = tau(se.C.2.SIZE.Y, n.C.2.SIZE.Y)

# clutch size (#eggs) per second nest (ADULTS)
mu.C.2.SIZE.A = 10.4		
se.C.2.SIZE.A = 0.71
n.C.2.SIZE.A  = 7

tau.C.2.SIZE.A = tau(se.C.2.SIZE.A, n.C.2.SIZE.A)

# nest survival probability for first nest (YEARLINGS)
mu.N.SURV.1.Y = 0.16 	
se.N.SURV.1.Y = 0.05
n.N.SURV.1.Y  = 37

tau.N.SURV.1.Y = tau(se.N.SURV.1.Y, n.N.SURV.1.Y)

# nest survival probability for first nest (ADULTS)
mu.N.SURV.1.A = 0.18 	
se.N.SURV.1.A = 0.05
n.N.SURV.1.A  = 35

tau.N.SURV.1.A = tau(se.N.SURV.1.A, n.N.SURV.1.A)

# nest survival probability for second nest (YEARLINGS)
mu.N.SURV.2.Y = 0.28 	
se.N.SURV.2.Y = 0.07
n.N.SURV.2.Y  = 12

tau.N.SURV.2.Y = tau(se.N.SURV.2.Y, n.N.SURV.2.Y)

# nest survival probability for second nest (ADULTS)
mu.N.SURV.2.A = 0.31 	
se.N.SURV.2.A = 0.08
n.N.SURV.2.A  = 8

tau.N.SURV.2.A = tau(se.N.SURV.2.A, n.N.SURV.2.A)

# hatching success (Chicks / Eggs)
mu.HS         = 0.80			
se.HS         = 0.05
n.HS          = 28

tau.HS = tau(se.HS, n.HS)

# brood survival probability (YEARLINGS)
mu.B.SURV.Y   = 0.34		
se.B.SURV.Y   = 0.07
n.B.SURV.Y    = 15

tau.B.SURV.Y = tau(se.B.SURV.Y, n.B.SURV.Y)

# brood survival probability (ADULTS)
mu.B.SURV.A   = 0.34		
se.B.SURV.A   = 0.07
n.B.SURV.A    = 20

tau.B.SURV.A = tau(se.B.SURV.A, n.B.SURV.A)

# fledging success (Fledglings / Chicks)
mu.FS         = 0.48			
se.FS         = 0.04
n.FS          = 16

tau.FS = tau(se.FS, n.FS)

# Juvenile Female survival
mu.S.J        = 0.38
se.S.J        = 0.002
n.S.J         = 18

tau.S.J = tau(se.S.J, n.S.J)

# Yearling Female survival
mu.S.Y        = 0.34			
se.S.Y        = 0.001
n.S.Y         = 53

tau.S.Y = tau(se.S.Y, n.S.Y)

# Adult Female survival
mu.S.A        = 0.42			
se.S.A        = 0.002
n.S.A         = 63

tau.S.A = tau(se.S.A, n.S.A)

((N.PROB * mu.C.1.SIZE.Y * mu.N.SURV.1.Y) +
            (1 - mu.N.SURV.1.Y) * (RE.N.PROB * mu.C.2.SIZE.Y * mu.N.SURV.2.Y)) *
          mu.HS * mu.B.SURV.Y * mu.FS * 0.5

# # Fecundity of yearlings
# F.Y = ((N.Prob * C.1.Size * N.Surv.1.Y) +
# 		(1 - N.Surv.1.Y) *
# 		(Re.N.Prob * C.2.Size * N.Surv.2.Y)) *
# 	  HS * B.Surv * FS * 0.5

# # Fecundity of adults
# F.A = ((N.Prob * C.1.Size * N.Surv.1.A) +
# 		 (1 - N.Surv.1.A) *
# 		 (Re.N.Prob * C.2.Size * N.Surv.2.A)) *
# 	  HS * B.Surv * FS * 0.5

# # Translate to matrix to determine eigenvalues->lambda
# values.MAT <- list(F.Y = F.Y,
# 					F.A = F.A,
# 					S.J = S.J,
# 					S.Y = S.Y,
# 					S.A = S.A)

# elements.MAT <- expression(F.Y * S.J, F.A * S.J,
# 						    S.Y   ,    S.A)

# Vrs.MAT <- sapply(elements.MAT, eval, values.MAT, NULL)

# A.MAT <- matrix(Vrs.MAT, nrow=2, byrow=TRUE)

# evs <- eigen(A.MAT)

# # LAMBDA
# imax <- which.max(evs$values)
# LAMBDA.VR <- Re(evs$values[imax])

# --------------------------------
# Initialize COUNTS
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

leks$date <- year(mdy(leks$date))
leks$n.tot[leks$n.tot=="#VALUE!"] <- NA
leks$n.m[leks$n.m=="#VALUE!"] <- NA
leks$n.f[leks$n.f=="#VALUE!"] <- NA

leks$n.tot <- as.numeric(leks$n.tot)
leks$n.m <- as.numeric(leks$n.m)
leks$n.f <- as.numeric(leks$n.f)

# Lek count data (broadly following Winder et al. 2015)
# (1) Find discount ratio for non sex-ratio counts (includes both trap and flush counts) 
#		from sex-ratio counts (includes both trap and flush counts)
# (2) Discount all the counts (trap and flush) that don't have sex-differentiation
# (3) Find the NUMBER OF COUNTS of each type for each LEK and YEAR
# (4) Find the maximum count (post-discount) of each type for each year
# (5) Calculated a weighted mean for max trap and flush counts by the proportion of each count type
# (6) Convert to female counts [currently n.m = n.f]
# (7) Assemble female count matrix
# (8) Calculate population size across all leks, and proportion of leks surveyed

# (1)
sex_ratio_discount <- leks %>%
						subset(!is.na(n.f) & !is.na(n.m)) %>%			# take the counts that have m/f differentiation
						transmute(coeff = n.m / n.tot) %>%				# calculate the male representation
						summarise(mean = mean(coeff, na.rm=TRUE), 		# remove zero counts (which give a NaN ratio)
								  sd = sd(coeff, na.rm=TRUE))			# summarize coefficient values

discount <- function(n.tot, n.m) {
	adjustedCount <- n.m
	coeff <- sex_ratio_discount$mean

	if (is.na(n.m)) {								# if n.m is NOT N/A (i.e., there was a sex count), just keep the count
		adjustedCount <- (n.tot * coeff)			
	} 

	return(adjustedCount)
}

# (2)
adjusted_male_counts <- leks %>%
							select(trapped, id, date, n.tot, n.m) %>% 
							mutate(adjusted.n.m = pmap_dbl(list(n.tot, n.m), discount))	# discounts at every row


# (3)
count_numbers <- adjusted_male_counts %>%
					group_by(id, date, trapped) %>%							# for each lek and year
					dplyr::count() %>%										# how many counts of each type?
					spread(trapped, n, fill=0) %>%									# each lek and year with just a single column
					set_colnames(c("id", "date", "n.flush", "n.trap")) %>%	# rename columns
					mutate(total.counts = n.flush+n.trap) %>%				# convert to proportions
					mutate(p.flush=n.flush/total.counts, 
						   p.trap=n.trap/total.counts) %>%
					select(id, date, p.flush, p.trap)						# extract just lek, date, proportions

# (4)+(5) NOTE ROUNDING!!
weightedMean <- function(n.fl, n.tr, p.fl, p.tr) {
	if (is.na(n.fl)) { n.fl <- 0 }		# this removes the NAs but doesn't matter because p will equal 0 for that count
	if (is.na(n.tr)) { n.tr <- 0 }

	weighted.n <- (n.fl * p.fl) + (n.tr * p.tr)
	return(round(weighted.n))						# note rounding!!!
}

max_counts <- adjusted_male_counts %>%
					group_by(id, date, trapped) %>%				# for each lek, year, and count type
					top_n(1, adjusted.n.m) %>%					# what's the max adjusted male count?
					sample_n(1) %>%								# top_n returns ties with no option (silly imo), so just pick one 
					select(id, date, trapped, n=adjusted.n.m) %>%
					spread(trapped, n) %>%
					set_colnames(c("id", "date", "n.m.flush", "n.m.trap")) %>%
					ungroup() %>%
					left_join(count_numbers, by=c("id","date")) %>%				# join with the trap proportion values
					mutate(final.n.m = pmap_dbl(list(n.m.flush, n.m.trap,
													 p.flush, p.trap), 
												weightedMean))					# calculate final weighted mean NOTE ROUNDED!

# (6) NOTE AGAIN N.M = N.F, no further adjustements right now

# (7) Convert df to matrix
female.counts <- max_counts %>%
					select(id, date, final.n.m) %>%
					spread(date, final.n.m) %>%
					select(-id) %>%
					as.matrix()

# (8)
N.LEKS <- nrow(female.counts)
N.YEARS <- ncol(female.counts)

P.SURVEYED <- rep(0, times=N.YEARS)
for (t in 1:N.YEARS) {
	num.surveyed <- length(which(!is.na(female.counts[,t])))
	P.SURVEYED[t] <- num.surveyed / N.LEKS
}

COUNTS <- colSums(female.counts, na.rm=TRUE)