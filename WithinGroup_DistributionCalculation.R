
# ============================== %
# Calculate distribution position of score within group %
# Author: Stephanie Kirmer
# July 2016
# This script was really fun to write - it takes a score, compares it to the scores of the rest of the group, and produces the distribution point for the score.
# There was some added complexity because the groups were identified by prior year.
# ============================== %

# Load libraries ####
library("foreign")
library("sas7bdat")
library(dplyr)
library(data.table)
library(reshape2)
library("readstata13")
library(sem)
library(doBy)
library(ggplot2)
require(sem)
require(AER)
require(ivpack)
require(lubridate)
library(mondate)
require("ggrepel")



## ===============      WITHIN SCHOOL DISTRIBUTION     ================ ####

#=====================================##=====================================#

test1a<-#load your data
test2a<-#load your data
test3a<-#load your data
test4a<-#load your data

master2013 <- read.dta13("/master2013.dta")
master2014 <- read.dta13("/master2014.dta")
master2015 <- read.dta13("/master2015.dta")

# Trim down to needed pieces (pulling together the district-wide datasets for test scores)
test1 <- test1a[test1a$Discipline=="Mathematics" & test1a$TermName=="Spring 2013-2014" & test1a$Grade==8,
                  c("TermName","sid","SchoolID","SchoolName","StudentEthnicGroup","StudentGender","Grade",
                    "Discipline","TestType","TestName","TestPercentile","SchoolYear")]
test3 <- test3a[test3a$grade==9,c("sid","yr","grade","gender","race","disab","lunch","mathxil","compxil")]
test4 <- test4a[test4a$grade==8,c("sid","yr","grade","gender","race","disab","lunch","mathxil","compxil", "schlid")]
test2 <- test2a[test2a$grade=="09",c("sid","yr","grade","schoolid","gender","race","lunch","mathxil","compxil","school_year")]

#Minor data cleaning
test3x <- merge(test3, master2013[,c("sid", "schlid")], by.x="sid", by.y="sid", all.x=T)
test2$schlid <- test2$schoolid
test1 <-  rename(test1, schlid=SchoolID)
test1$mathxil <- as.numeric(test1$TestPercentile)

#Just putting together the student/school records- attaching the current year's school with the prior year's test score.
#We are comparing the kids to their current peers but on the PRIOR year (baseline) test score, so this is IMPORTANT.
test4.s <- merge(test4, master2014[,c("sid", "schlid")], by.x="sid", by.y="sid")
test3.s <- merge(test3x, master2014[,c("sid", "schlid")], by.x="sid", by.y="sid")
test2.s <- merge(test2, master2015[,c("sid", "schlid")], by.x="sid", by.y="sid")
test1.s <- merge(test1, master2015[,c("sid", "schlid")], by.x="sid", by.y="sid")

test4.s <-  rename(test4.s, school_post1=schlid.y, school_pre=schlid.x)
test3.s <-  rename(test3.s, school_post1=schlid.y, school_pre=schlid.x)
test2.s <-  rename(test2.s, school_post1=schlid.y, school_pre=schlid.x)
test1.s <-  rename(test1.s, school_post1=schlid.y, school_pre=schlid.x)

# Make a table of all the test/subject/year options for placing in the function below. ####
datatest1 <- expand.grid(testfile=c("test4.s"),masterfile="master_dataset",testcolumn=c("test4_pre"),scorecolumn=c("mathxil_pre"),yr=c(1))
datatest2 <- expand.grid(testfile=c("test3.s"),masterfile="master_dataset",testcolumn=c("test3_pre"),scorecolumn=c("mathxil_pre"),yr=c(1))
datatest3 <- expand.grid(testfile=c("test1.s"),masterfile="master_dataset",testcolumn=c("test1_pre"),scorecolumn=c("mathxil_pre"),yr=c(2))
datatest4 <- expand.grid(testfile=c("test2.s"),masterfile="master_dataset",testcolumn=c("test3_pre"),scorecolumn=c("mathxil_pre"),yr=c(2))

datatest <- rbind(datatest1, datatest2, datatest3, datatest4)

#quick formatting stuff
datatest$testfile <- as.character(datatest$testfile)
datatest$masterfile <- as.character(datatest$masterfile)
datatest$testcolumn <- as.character(datatest$testcolumn)

# Function to create the within-school distributions ####
within_school_score_dist <- function(testfile, masterfile, testcolumn, scorecolumn, schlidcolumn, yr){

  d <- list() #Create container for output data
  masterfile2 <- eval(parse(text=masterfile)) #Make the name of the student file recognized as an object
  testfile2 <- eval(parse(text=testfile)) #Make the name of the all-district scores file recognized as an object
  sidlist <- masterfile2[,"sid"]
  post_school <- testfile2[,"school_post1"]
  mathscore <- testfile2[,"mathxil"]

  for (i in sidlist){
    rownum <- which(sidlist==i) #Find out the student's row
    kidschool <- as.numeric(masterfile2[rownum,schlidcolumn]) #Find out which school we are looking at
    kidstudy <- masterfile2[rownum,"study"] #Find out which study the kid is in

    filtered_testfile <- testfile2[c(post_school==kidschool & !is.na(mathscore)),]
    #Filter for the correct school and correct test subject (we already know it's the right test)

    comparison_grp <- nrow(filtered_testfile) # How many kids in our comparison set?
    kidscore <- as.numeric(masterfile2[rownum,scorecolumn]) #Identify the student's score
    kidtest <- masterfile2[rownum, testcolumn] #Identify which test we are dealing with

    if (comparison_grp!=0){ #Make sure the comparison group is nonzero
      if (!(is.na(kidscore))) { #Make sure our kid has a score
        if (!(is.na(kidtest))) { #Make sure the kid took the test we are expecting
          if(kidstudy == yr){ #Student's study year matches the inputs (don't look for the wrong study for the year the kid was in)
            quan<-ecdf(filtered_testfile$mathxil) #Run the distribution within school for the math xil for the test.
            xil <- quan(kidscore) #Place the student's score within the school distro
            data <- data.frame(Row=rownum, SID=i, Score=kidscore, Percentile=xil, School=kidschool, Test=testcolumn, Study=yr)
            #Add what you want to see in your file to an object
            print(data)
            d[[i]] <- data #Place that object in your container list created above
          }
        }
      }
    }
  }
  assign(paste0("withinschoolscores_", testcolumn,"_", scorecolumn, "_", yr), envir=.GlobalEnv, rbindlist(d)) #Create the output with dynamic naming
}
mapply(within_school_score_dist, datatest$testfile, datatest$masterfile, datatest$testcolumn, datatest$scorecolumn, "schlid", datatest$yr) #Run the function

# Merge the output files together ####

withinschool_pctile_pre <- rbind(withinschoolscores_test4_pre_mathxil_pre_1[,c("SID", "Percentile", "School"), with=FALSE],
                                 withinschoolscores_test3_pre_mathxil_pre_1[,c("SID", "Percentile", "School"), with=FALSE],
                                 withinschoolscores_test3_pre_mathxil_pre_2[,c("SID", "Percentile", "School"), with=FALSE],
                                 withinschoolscores_test1_pre_mathxil_pre_2[,c("SID", "Percentile", "School"), with=FALSE])

withinschool_pctile_pre <-  plyr::rename(withinschool_pctile_pre,c("Percentile"="withinschool_math_dist_xil", "School"="test_dist_schlid"))

#Pull everything into the original dataset
master_dataset <- merge(master_dataset, withinschool_pctile_pre, by.x="sid", by.y="SID", all.x=T)
