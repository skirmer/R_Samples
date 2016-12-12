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
sd(dataset$test, na.rm=T)
sd(dataset$gpa, na.rm=T)

#Create the different sub-populations- choose the one you need in a given moment, don't run both or you'll just get the last one.
population <- filter(dataset, test > (20-4.5), test < (21+4.5), gpa > (3-.75), gpa < (3+.75))

#trait1
t1 <- prop.table(table(population$trait1,population$trait2, useNA="ifany", dnn=c("t1", "t2")),2)
t1_n <- table(population$trait2, useNA="ifany")
t1_t2 <- rbind(t1[2,], t1_n)

#Look at the results
rownames(t1_t2) <- c("t1_prob", "n")

# Run the power analysis, inputting different variations of the results - you can use the function help section to decide what you want to do there.
power.prop.test(n=150, p1=.264, p2=NULL, sig.level=.05, power=.8)

