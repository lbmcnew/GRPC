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
# Mean pre- and post-construction vital rates
# FROM Table 4.3 Sandercock DE FINAL TECHNICAL REPORT 
# --------------------------------

# given a mean vital rate value for a pre and post period,
# return a timeseries for the given three-year periods
# Pre-construction = 2007, 2008
# Post-construction = 2009, 2010, 2011
tsVRS <- function(pre, post) {
	ts <- c(rep(pre, times=2), rep(post, times=3))
	return(ts)
}

# NOTE set to 1 following McNew et al. 2012 methods
N.PROB        = 1 	# nesting probability per female
RE.N.PROB     = 1	# re-nesting probability per female

# clutch size (#eggs) per first nest
pre.C.1.SIZE  = 12.3
post.C.1.SIZE = 12.8
C.1.SIZE = tsVRS(pre.C.1.SIZE, post.C.1.SIZE)

# clutch size (#eggs) per second nest
pre.C.2.SIZE  = 10.7
post.C.2.SIZE = 10.6
C.2.SIZE = tsVRS(pre.C.2.SIZE, post.C.2.SIZE)

# nest survival probability for first nest
pre.N.SURV.1  = 0.25
post.N.SURV.1 = 0.31
N.SURV.1 = tsVRS(pre.N.SURV.1, post.N.SURV.1)

# nest survival probability for second nest
# TODO nest survival vs. 35 day nest success in Table 4.3?
pre.N.SURV.2  = 0.36
post.N.SURV.2 = 0.30
N.SURV.2 = tsVRS(pre.N.SURV.2, post.N.SURV.2)

# hatching success for first nests (Chicks / Eggs)
pre.HS.1  = 0.81
post.HS.1 = 0.87
HS.1 = tsVRS(pre.HS.1, post.HS.1)

# hatching success for second nests (Chicks / Eggs)
pre.HS.2  = 0.81
post.HS.2 = 0.78
HS.2 = tsVRS(pre.HS.2, post.HS.2)

# brood survival probability (all 25 days)
pre.B.SURV  = 0.50
post.B.SURV = 0.38
B.SURV = tsVRS(pre.B.SURV, post.B.SURV)

# fledging success (Fledglings / Chicks; all 25 days)
pre.FS  = 0.58
post.FS = 0.35
FS = tsVRS(pre.FS, post.FS)

# Juvenile Female survival
# NOTE cannot find this value in the DoE report,
# 	So using McNew et al. 2012 [pre] values for both pre- and post-
pre.S.J = 0.38
post.S.J = 0.38
S.J = tsVRS(pre.S.J, post.S.J)

# Adult Female survival
# FROM cumulative mortality reports (Sandercock 2013, pg. 51)
pre.S.A  = 0.274
post.S.A = 0.543
S.A = tsVRS(pre.S.A, post.S.A) 

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
			subset(!is.na(n.f) & !is.na(n.m)) %>%	   # take the counts that have m/f differentiation
			transmute(coeff = n.m / n.tot) %>%	   # calculate the male representation
			summarise(mean = mean(coeff, na.rm=TRUE),  # remove zero counts (which give a NaN ratio)
			sd = sd(coeff, na.rm=TRUE))		   # summarize coefficient values

discount <- function(n.tot, n.m) 
{
	adjustedCount <- n.m
	coeff <- sex_ratio_discount$mean

	if (is.na(n.m)) {			  # if n.m is NOT N/A (i.e., there was a sex count), just keep the count
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
		 	group_by(id, date, trapped) %>%				# for each lek and year
		 	dplyr::count() %>%					# how many counts of each type?
		 	spread(trapped, n, fill=0) %>%			 	# each lek and year with just a single column
		 	set_colnames(c("id", "date", "n.flush", "n.trap")) %>%	# rename columns
		 	mutate(total.counts = n.flush+n.trap) %>%		# convert to proportions
		 	mutate(p.flush=n.flush/total.counts, 
			   p.trap=n.trap/total.counts) %>%
		 	select(id, date, p.flush, p.trap)			# extract just lek, date, proportions

# (4)+(5) NOTE ROUNDING!!
weightedMean <- function(n.fl, n.tr, p.fl, p.tr) 
{
	if (is.na(n.fl)) { n.fl <- 0 }		# this removes the NAs but doesn't matter because p will equal 0 for that count
	if (is.na(n.tr)) { n.tr <- 0 }

	weighted.n <- (n.fl * p.fl) + (n.tr * p.tr)
	return(round(weighted.n))						# note rounding!!!
}

max_counts <- adjusted_male_counts %>%
	      	group_by(id, date, trapped) %>%				 	# for each lek, year, and count type
	      	top_n(1, adjusted.n.m) %>%					# what's the max adjusted male count?
	      	sample_n(1) %>%							# top_n returns ties with no option (silly imo), so just pick one 
	      	select(id, date, trapped, n=adjusted.n.m) %>%
	      	spread(trapped, n) %>%
	      	set_colnames(c("id", "date", "n.m.flush", "n.m.trap")) %>%
	      	ungroup() %>%
	      	left_join(count_numbers, by=c("id","date")) %>%			# join with the trap proportion values
	      	mutate(final.n.m = pmap_dbl(list(n.m.flush, n.m.trap,
						 p.flush, p.trap), 
						 weightedMean))			# calculate final weighted mean NOTE ROUNDED!

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