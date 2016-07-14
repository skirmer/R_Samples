# Power Analysis template #
# Stephanie Kirmer #
# July 2016 #
# These steps help you learn what N you need to detect a certain effect size in an RCT before you proceed. #

library("stats")
library(mondate)
library(data.table)
library(dplyr)
library(reshape2)
library(sas7bdat)
library(rgeos)
library(rgdal)
library(ggplot2)
library(sp)
library(spdep)
library(raster)
library(RColorBrewer)
library(maptools)
library(classInt)
library(lubridate)
library(reshape)

#Pull in whatever dataset you need. This uses a group of students and examines some educational outcomes.

#Quick review:
#power: pct of attempts that result in a correct rejection of the null hypothesis
#p1 vs p2: difference is the size of effect you're detecting
#n: size of sample you need to get to detect the effect shown at the frequency your power indicates
#sig.level: the statistical significance you can expect for your estimate


#Create the outcomes if you haven't yet - looking for a binary variable here.

#check the SD to make sure thresholds are right for your subpopulations below
sd(dataset$act, na.rm=T)
sd(dataset$gpa_12th, na.rm=T)

#Create the different sub-populations- choose the one you need in a given moment, don't run both or you'll just get the last one.
population <- filter(dataset, act > (20-4.5), act < (21+4.5), gpa_12th > (3-.75), gpa_12th < (3+.75),  toupper(STUDENT_RACE)=="HISPANIC")
population <- filter(dataset, act > (20-4.5), act < (20+4.5), gpa_12th > (3-.75), gpa_12th < (3+.75))

#anycollege
anycol <- prop.table(table(population$went_college_flag,population$STUDENT_FOODSERVICE_INDICATOR, useNA="ifany", dnn=c("college", "frl")),2)
anycol_n <- table(population$STUDENT_FOODSERVICE_INDICATOR, useNA="ifany")
anycol_frl <- rbind(anycol[2,], anycol_n)

anycol <- prop.table(table(population$went_college_flag,population$gender, useNA="ifany", dnn=c("college", "gend")),2)
anycol_n <- table(population$opi_gender, useNA="ifany")
anycol_gend <- rbind(anycol[2,], anycol_n)

#graduation
population <- filter(dataset, act > (20-4.5), act < (21+4.5), gpa_12th > (3-.75), gpa_12th < (3+.75),  toupper(STUDENT_RACE)=="HISPANIC", cohort=="2008")
population <- filter(dataset, act > (20-4.5), act < (20+4.5), cohort=="2008", gpa_12th > (3-.75), gpa_12th < (3+.75))

graduation <- prop.table(table(population$GRADUATED,population$STUDENT_FOODSERVICE_INDICATOR, useNA="ifany", dnn=c("grad", "frl")),2)
graduation_n <- table(population$STUDENT_FOODSERVICE_INDICATOR, useNA="ifany")
graduation_frl <- rbind(graduation[2,], graduation_n)

graduation <- prop.table(table(population$GRADUATED,population$gender, useNA="ifany", dnn=c("grad", "gend")),2)
graduation_n <- table(population$gender, useNA="ifany")
graduation_gend <- rbind(graduation[2,], graduation_n)


#Look at the results
rownames(anycol_frl) <- c("coll_prob", "n")
rownames(graduation_frl) <- c("grad_prob", "gr_n")
rbind(anycol_frl,graduation_frl)

rownames(anycol_gend) <- c("coll_prob", "n")
rownames(graduation_gend) <- c("grad_prob", "gr_n")
rbind(anycol_gend, graduation_gend)


# Run the power analysis, inputting different variations of the results - you can use the function help section to decide what you want to do there.
power.prop.test(n=150, p1=.264, p2=NULL, sig.level=.05, power=.8)

