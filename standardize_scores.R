# Standardize test scores to control mean#
# Stephanie Kirmer #
# July 2016 #

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
library(sem)
library(AER)
library(ivpack)
library(lubridate)
library(mondate)
library("ggrepel")


# Calculating scaled scores ####

#Make a table of all the test/subject/year options for placing in the function below.

test_scores_1 <- expand.grid(masterfile="master_dataset",
                             testcolumn=c("plan_pre", "explore_gr9_pre", "explore_gr8_pre", "nwea_pre"),
                             scorecolumn=c("readxil_pre","mathxil_pre"),
                             yr=c(1,2))

test_scores_2 <- expand.grid(masterfile="master_dataset",
                             testcolumn=c("plan_post1", "explore_post1", "act_post1"),
                             scorecolumn=c("readxil_post1","mathxil_post1"),
                             yr=c(1,2))

test_scores_3 <- expand.grid(masterfile="master_dataset",
                             testcolumn=c("plan_post2","explore_post2","act_post2"),
                             scorecolumn=c("readxil_post2","mathxil_post2"),
                             yr=c(1,2))

test_scores <- rbind(test_scores_1, test_scores_2, test_scores_3)

#tidy up
rm(test_scores1, test_scores2, test_scores3)

test_scores$masterfile <- as.character(test_scores$masterfile)
test_scores$testcolumn <- as.character(test_scores$testcolumn)
test_scores$scorecolumn <- as.character(test_scores$scorecolumn)

#This function takes test data and scales it to the control mean by test, subject, & year
scaled_scores <- function(masterfile, testcolumn, scorecolumn, yr){
  #Source data file, the name of the test column, and the name of the score column, and study number

  d <- list() #Create container for output data
  masterfile2 <- eval(parse(text=masterfile)) #Make the function recognize the name of the file as an object
  sidlist <- masterfile2[,"sid"]

  newtestfile <- masterfile2[c(masterfile2[,"dmatch"]==0
                               & masterfile2[,testcolumn]==1
                               & masterfile2[,"study"]==yr),]# Filter the file that will make the mean and sd #
  mean_score <- mean(newtestfile[,scorecolumn], na.rm=T) #Get the mean#
  sd_score <- sd(newtestfile[,scorecolumn], na.rm=T) #Get the sd#

  for (i in sidlist){
    rownumber <- which(sidlist==i) # identify the row #
    x <- as.numeric(masterfile2[rownumber,scorecolumn]) # get the score out #
    testval <- masterfile2[rownumber,testcolumn] #Find out which test it was

    try( #Using try to bypass any irrelevant warning messages
      if(!is.na(testval) & length(testval) > 0 & !is.na(x) & !is.na(mean_score)){ # If the test is correct, there's a mean to compare to, and a score is present, do the calculations#
        scaledversion <- (x-mean_score)/sd_score # Calculate the score
        data <- data.frame(Row=rownumber, SID=i, rawscore=x, Scaledscore=scaledversion, Test=testcolumn,
                           Scoretype=scorecolumn, Mean=mean_score, StDev=sd_score, Study=yr)
        print(data) #output the entire row of data that is produced for checking
        d[[i]] <- data #Add your row of data to list d
      }
    )
  }

  assign(paste0("outputfile_", testcolumn,"_", scorecolumn,"_", yr), envir=.GlobalEnv, rbindlist(d))
  #make the output file for this loop- paste together the name so it is specific to the inputs,
  #make sure it is created in the global environment, and bind together everything you have saved into list d to create it.
}

compiled_ss <- cmpfun(scaled_scores)
mapply(compiled_ss, test_scores$masterfile, test_scores$testcolumn, test_scores$scorecolumn, test_scores$yr) #Run the function

# Combine the various output frames from the function ####
scaled_scores_pre_read <- rbind(outputfile_plan_pre_readxil_pre_1[,c("SID", "Scaledscore"), with=FALSE],
                                outputfile_explore_gr8_pre_readxil_pre_1[,c("SID", "Scaledscore"), with=FALSE],
                                outputfile_explore_gr9_pre_readxil_pre_1[,c("SID", "Scaledscore"), with=FALSE],
                                outputfile_nwea_pre_readxil_pre[,c("SID", "Scaledscore"), with=FALSE])

scaled_scores_pre_read <- plyr::rename(scaled_scores_pre_read,c("Scaledscore"="readxil_z_pre"))

scaled_scores_pre_math <- rbind(outputfile_plan_pre_mathxil_pre_1[,c("SID", "Scaledscore"), with=FALSE],
                                outputfile_explore_gr8_pre_mathxil_pre_1[,c("SID", "Scaledscore"), with=FALSE],
                                outputfile_explore_gr9_pre_mathxil_pre_1[,c("SID", "Scaledscore"), with=FALSE],
                                outputfile_nwea_pre_mathxil_pre[,c("SID", "Scaledscore"), with=FALSE])

scaled_scores_pre_math <- plyr::rename(scaled_scores_pre_math,c("Scaledscore"="mathxil_z_pre"))


scaled_scores_post1_read <- rbind(outputfile_plan_post1_readxil_post1_1[,c("SID", "Scaledscore"), with=FALSE],
                                  outputfile_explore_post1_readxil_post1_1[,c("SID", "Scaledscore"), with=FALSE],
                                  outputfile_act_post1_readxil_post1[,c("SID", "Scaledscore"), with=FALSE])

scaled_scores_post1_read <- plyr::rename(scaled_scores_post1_read,c("Scaledscore"="readxil_z_post1"))

scaled_scores_post1_math <- rbind(outputfile_plan_post1_mathxil_post1_1[,c("SID", "Scaledscore"), with=FALSE],
                                  outputfile_explore_post1_mathxil_post1_1[,c("SID", "Scaledscore"), with=FALSE],
                                  outputfile_act_post1_mathxil_post1[,c("SID", "Scaledscore"), with=FALSE])

scaled_scores_post1_math <- plyr::rename(scaled_scores_post1_math,c("Scaledscore"="mathxil_z_post1"))


scaled_scores_post2_read <- rbind(
  outputfile_explore_post2_readxil_post2_1[,c("SID", "Scaledscore"), with=FALSE],
  outputfile_act_post2_readxil_post2_1[,c("SID", "Scaledscore"), with=FALSE],
  outputfile_plan_post2_readxil_post2_1[,c("SID", "Scaledscore"), with=FALSE],
scaled_scores_post2_read <- plyr::rename(scaled_scores_post2_read,c("Scaledscore"="readxil_z_post2"))

scaled_scores_post2_math <- rbind(
  outputfile_explore_post2_mathxil_post2_1[,c("SID", "Scaledscore"), with=FALSE],
  outputfile_act_post2_mathxil_post2_1[,c("SID", "Scaledscore"), with=FALSE],
  outputfile_plan_post2_mathxil_post2_1[,c("SID", "Scaledscore"), with=FALSE],

scaled_scores_post2_math <- plyr::rename(scaled_scores_post2_math,c("Scaledscore"="mathxil_z_post2"))

scaled_scores_complete <- merge(scaled_scores_pre_read, scaled_scores_pre_math, by="SID", all.x=T, all.y=T)
scaled_scores_complete <- merge(scaled_scores_complete, scaled_scores_post1_read, by="SID", all.x=T, all.y=T)
scaled_scores_complete <- merge(scaled_scores_complete, scaled_scores_post1_math, by="SID", all.x=T, all.y=T)
scaled_scores_complete <- merge(scaled_scores_complete, scaled_scores_post2_read, by="SID", all.x=T, all.y=T)
scaled_scores_complete <- merge(scaled_scores_complete, scaled_scores_post2_math, by="SID", all.x=T, all.y=T)

#Pull everything into the original dataset
master_dataset <- merge(master_dataset, scaled_scores_complete, by.x="sid", by.y="SID", all.x=T)