# setwd("~/GRPC")
setwd("C://Users//Liam//Documents//GRPC") 

require(tidyverse)
require(lubridate)
require(magrittr)
require(jagsUI)

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
LAMBDA.VR <- Re(evs$values[imax])


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
					sample_n(1) %>%								# top_n returns ties with no option (#DUMB), so just pick one 
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