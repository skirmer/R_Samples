# Permutation Tests (Sensitivity Analysis) Procedure #
# Stephanie Kirmer #
# August 2016 #

#This project involved running a 2-stage least squares regression estimating the effect of an academic intervention
#on a number of outcomes. In order to check the sensitivity of the analysis, we used permutation tests to re-sample our
#study population, abiding by the original randomization block and control/treatment divide sizes.
#Literally, we are taking our same data and reshuffling it 100,000 times to see if the results would have been
#different if our treatment and randomization assignments were different.

#A good resource for more information on the basics: http://thomasleeper.com/Rcourse/Tutorials/permutationtests.html


# Steps: 100,000 resamplings following the original sampling design
# Sample with replacement, keeping the blocks balanced
# Recalculate the standardized test scores pre and post, reading and math
# Calculate regression model, keep t-stat for the ITT on the main outcomes (or p-value if you prefer)
# Save out to file


#Part of what is unique about this project is that we used control-mean standardization on test scores throughout. 
#This means that each permutation must re-standardize the scores for pre-tests, and when post-test is the outcome, that has to be done too.
#There are two functions here, one that does post-test and one that doesn't, in order to save some time in the code running.
#It would also be an option to make this all a singular function and use conditional language to choose to include the post-test standardization.

#This code also calculates the time of the system processing, in case you're interested in that.
#dmatch represents the treatment/control assignment, and blocknum represents the randomization block.

library(resample)
library(boot)
library(data.table)

master_dataset_1<- read.csv("/export/projects/BAM/Match_AnalysisFiles_May2016/SourceFiles/Archive/analysisdata_study1.csv", stringsAsFactors=F)
# print(table(master_dataset_1$blocknum, master_dataset_1$dmatch ,dnn=c("block", "pseudoT"), useNA="ifany"))


# function to obtain regression output for each resample ####

bs_nopost <- function(data, indices, M) {

  d <- data[indices,] # allows boot to select sample
  # Create variables to identify size of block and size of block/assignment
  d <- as.data.table(d, keep.rownames=FALSE)
  blocktxsizes<- summarize(group_by(d, blocknum),
                           size = sum(dmatch==1, na.rm=T))
  d <- merge(d, blocktxsizes, by.x="blocknum", by.y="blocknum", all.x=T)
#   table(d$size, d$blocknum) #Make sure block/tx sizes are right

  #Sort students within schools by this number and select first x (number of treated in school)

  #Identify the "pseudoT" group
  set.seed(145560) #For ensuring replicability
  d$randomnum <- sample(1:nrow(d))
  #print(summary(d$randomnum[d$blocknum==20])) #Check to make sure the random number assignment is replicable
  d[, blkrank:= rank(-randomnum), by=blocknum]
  d[blkrank <= size, pseudoT:=1]
  d[blkrank > size, pseudoT:=0]

  # Standardize the test scores (baseline and post 1) so that the model is right
  # Going with the simplistic tedious approach -  this is also a bit faster than the alternatives, but 
  mean_math_plan_pre <- mean(d$mathxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  sd_math_plan_pre <- sd(d$mathxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  mean_math_exp9_pre <- mean(d$mathxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  sd_math_exp9_pre <- sd(d$mathxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  mean_math_exp8_pre <- mean(d$mathxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)
  sd_math_exp8_pre <- sd(d$mathxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)

  mean_read_plan_pre <- mean(d$readxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  sd_read_plan_pre <- sd(d$readxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  mean_read_exp9_pre <- mean(d$readxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  sd_read_exp9_pre <- sd(d$readxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  mean_read_exp8_pre <- mean(d$readxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)
  sd_read_exp8_pre <- sd(d$readxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)

  #if anything is NA then make it zero - you'll get a message if this is an issue.
  sd_math_plan_pre[is.na(sd_math_plan_pre)] <- 0


  # Calculate new standardization #
  d[d$plan_pre==1, rescaled_math_plan_pre := (d$mathxil_pre-mean_math_plan_pre)/sd_math_plan_pre]
  d[d$explore_gr9_pre==1, rescaled_math_exp9_pre := (d$mathxil_pre-mean_math_exp9_pre)/sd_math_exp9_pre]
  d[d$explore_gr8_pre==1, rescaled_math_exp8_pre := (d$mathxil_pre-mean_math_exp8_pre)/sd_math_exp8_pre]

  d[d$plan_pre==1, rescaled_read_plan_pre := (d$readxil_pre-mean_read_plan_pre)/sd_read_plan_pre]
  d[d$explore_gr9_pre==1, rescaled_read_exp9_pre := (d$readxil_pre-mean_read_exp9_pre)/sd_read_exp9_pre]
  d[d$explore_gr8_pre==1, rescaled_read_exp8_pre := (d$readxil_pre-mean_read_exp8_pre)/sd_read_exp8_pre]

  d[,rescaled_math_pre := rescaled_math_plan_pre]
  d[is.na(rescaled_math_pre),rescaled_math_pre := rescaled_math_exp9_pre]
  d[is.na(rescaled_math_pre),rescaled_math_pre := rescaled_math_exp8_pre]
  d[is.na(rescaled_math_pre),rescaled_math_pre := 0]

  d[,rescaled_read_pre := rescaled_read_plan_pre]
  d[is.na(rescaled_read_pre),rescaled_read_pre := rescaled_read_exp9_pre]
  d[is.na(rescaled_read_pre),rescaled_read_pre := rescaled_read_exp8_pre]
  d[is.na(rescaled_read_pre),rescaled_read_pre := 0]

  d <- as.data.frame(d)

  # Run the actual LM and return the T-statistic for pseudoT(dmatch proxy)

I <- eval(parse(text=M))

  runmodel <- function(I) {
    fit <- lm(I~pseudoT+ blocknum+ d13andunder+d14+d15+d16+d17andover+ dlearningdisabled+dfreelunch+ dblack+ dhispanic+dother+
               dgrade9+ dgrade10+ gpa_pre_zeros+ numAs_pre+  numBs_pre+ numCs_pre+ numDs_pre+ numFs_pre+ missing_gpa_pre+ days_absent_pre_zeros+
               missing_attend_pre+ rescaled_math_pre+  rescaled_read_pre +mathxil_z_pre_missing+readxil_z_pre_missing+ oss_dis_pre_zeros+
               incidents_pre_zeros+any_arrests_pre+violent_pre+property_pre+drug_pre, data=d)
    summary(fit)
    #t-test
    output <<- summary(fit)$coefficients["pseudoT",3]
    #pval: #output <<- summary(fit)$coefficients["pseudoT",4]
  }
  runmodel(I)
  return(output)
}

bs_post <- function(data, indices, M) {

  d <- data[indices,] # allows boot to select sample
  # Create variables to identify size of block and size of block/assignment
  d <- as.data.table(d, keep.rownames=FALSE)
  blocktxsizes<- summarize(group_by(d, blocknum),
                           size = sum(dmatch==1, na.rm=T))
  d <- merge(d, blocktxsizes, by.x="blocknum", by.y="blocknum", all.x=T)

  #Sort students within schools by this number and select first x (number of treated in school)
  #Identify the "pseudoT" group
  set.seed(145560) #For ensuring replicability
  d$randomnum <- sample(1:nrow(d))

  d[, blkrank:= rank(-randomnum), by=blocknum]
  d[blkrank <= size, pseudoT:=1]
  d[blkrank > size, pseudoT:=0]

  # Standardize the test scores (baseline and post 1) so that the model is right
  # Going with the simplistic tedious approach
  mean_math_plan_pre <- mean(d$mathxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  sd_math_plan_pre <- sd(d$mathxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  mean_math_exp9_pre <- mean(d$mathxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  sd_math_exp9_pre <- sd(d$mathxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  mean_math_exp8_pre <- mean(d$mathxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)
  sd_math_exp8_pre <- sd(d$mathxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)

  mean_read_plan_pre <- mean(d$readxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  sd_read_plan_pre <- sd(d$readxil_pre[d$pseudoT==0 & d$plan_pre==1], na.rm=T)
  mean_read_exp9_pre <- mean(d$readxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  sd_read_exp9_pre <- sd(d$readxil_pre[d$pseudoT==0 & d$explore_gr9_pre==1], na.rm=T)
  mean_read_exp8_pre <- mean(d$readxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)
  sd_read_exp8_pre <- sd(d$readxil_pre[d$pseudoT==0 & d$explore_gr8_pre==1], na.rm=T)

  mean_math_plan_post1 <- mean(d$mathxil_post1[d$pseudoT==0 & d$plan_post1==1], na.rm=T)
  sd_math_plan_post1 <- sd(d$mathxil_post1[d$pseudoT==0 & d$plan_post1==1], na.rm=T)
  mean_math_exp_post1 <- mean(d$mathxil_post1[d$pseudoT==0 & d$explore_post1==1], na.rm=T)
  sd_math_exp_post1 <- sd(d$mathxil_post1[d$pseudoT==0 & d$explore_post1==1], na.rm=T)
  mean_math_act_post1 <- mean(d$mathxil_post1[d$pseudoT==0 & d$act_post1==1], na.rm=T)
  sd_math_act_post1 <- sd(d$mathxil_post1[d$pseudoT==0 & d$act_post1==1], na.rm=T)

  mean_read_plan_post1 <- mean(d$readxil_post1[d$pseudoT==0 & d$plan_post1==1], na.rm=T)
  sd_read_plan_post1 <- sd(d$readxil_post1[d$pseudoT==0 & d$plan_post1==1], na.rm=T)
  mean_read_exp_post1 <- mean(d$readxil_post1[d$pseudoT==0 & d$explore_post1==1], na.rm=T)
  sd_read_exp_post1 <- sd(d$readxil_post1[d$pseudoT==0 & d$explore_post1==1], na.rm=T)
  mean_read_act_post1 <- mean(d$readxil_post1[d$pseudoT==0 & d$act_post1==1], na.rm=T)
  sd_read_act_post1 <- sd(d$readxil_post1[d$pseudoT==0 & d$act_post1==1], na.rm=T)

  #if anything is NA then make it zero -  act has been a problem on this, might also run into it on plan pre
  sd_math_plan_pre[is.na(sd_math_plan_pre)] <- 0
  sd_math_act_post1[is.na(sd_math_act_post1)] <- 0

  # Calculate new standardization #
  d[d$plan_pre==1, rescaled_math_plan_pre := (d$mathxil_pre-mean_math_plan_pre)/sd_math_plan_pre]
  d[d$explore_gr9_pre==1, rescaled_math_exp9_pre := (d$mathxil_pre-mean_math_exp9_pre)/sd_math_exp9_pre]
  d[d$explore_gr8_pre==1, rescaled_math_exp8_pre := (d$mathxil_pre-mean_math_exp8_pre)/sd_math_exp8_pre]

  d[d$plan_pre==1, rescaled_read_plan_pre := (d$readxil_pre-mean_read_plan_pre)/sd_read_plan_pre]
  d[d$explore_gr9_pre==1, rescaled_read_exp9_pre := (d$readxil_pre-mean_read_exp9_pre)/sd_read_exp9_pre]
  d[d$explore_gr8_pre==1, rescaled_read_exp8_pre := (d$readxil_pre-mean_read_exp8_pre)/sd_read_exp8_pre]

  d[d$plan_post1==1, rescaled_math_plan_post1 := (d$mathxil_post1-mean_math_plan_post1)/sd_math_plan_post1]
  d[d$explore_post1==1, rescaled_math_exp_post1 := (d$mathxil_post1-mean_math_exp_post1)/sd_math_exp_post1]
  d[d$act_post1==1, rescaled_math_act_post1 := (d$mathxil_post1-mean_math_act_post1)/sd_math_act_post1]

  d[d$plan_post1==1, rescaled_read_plan_post1 := (d$readxil_post1-mean_read_plan_post1)/sd_read_plan_post1]
  d[d$explore_post1==1, rescaled_read_exp_post1 := (d$readxil_post1-mean_read_exp_post1)/sd_read_exp_post1]
  d[d$act_post1==1, rescaled_read_act_post1 := (d$readxil_post1-mean_read_act_post1)/sd_read_act_post1]

  d[,rescaled_math_pre := rescaled_math_plan_pre]
  d[is.na(rescaled_math_pre),rescaled_math_pre := rescaled_math_exp9_pre]
  d[is.na(rescaled_math_pre),rescaled_math_pre := rescaled_math_exp8_pre]
  d[is.na(rescaled_math_pre),rescaled_math_pre := 0]

  d[,rescaled_read_pre := rescaled_read_plan_pre]
  d[is.na(rescaled_read_pre),rescaled_read_pre := rescaled_read_exp9_pre]
  d[is.na(rescaled_read_pre),rescaled_read_pre := rescaled_read_exp8_pre]
  d[is.na(rescaled_read_pre),rescaled_read_pre := 0]

  d[,rescaled_math_post1 := rescaled_math_plan_post1]
  d[is.na(rescaled_math_post1),rescaled_math_post1 := rescaled_math_exp_post1]
  d[is.na(rescaled_math_post1) & sd_math_act_post1 >0,rescaled_math_post1 := rescaled_math_act_post1]
  d[is.na(rescaled_math_post1), rescaled_math_post1 := 0]

  d[,rescaled_read_post1 := rescaled_read_plan_post1]
  d[is.na(rescaled_read_post1),rescaled_read_post1 := rescaled_read_exp_post1]
  d[is.na(rescaled_read_post1) & sd_read_act_post1 >0,rescaled_read_post1 := rescaled_read_act_post1]
  d[is.na(rescaled_read_post1), rescaled_read_post1 := 0]

  d <- as.data.frame(d)

  # Run the actual LM and return the T-statistic for pseudoT(dmatch proxy)

  I <- eval(parse(text=M))

  runmodel <- function(I) {
    fit <- lm(I~pseudoT+ blocknum+ d13andunder+d14+d15+d16+d17andover+ dlearningdisabled+dfreelunch+ dblack+ dhispanic+dother+ 
                dgrade9+ dgrade10+ gpa_pre_zeros+ numAs_pre+  numBs_pre+ numCs_pre+ numDs_pre+ numFs_pre+ missing_gpa_pre+ days_absent_pre_zeros+
                missing_attend_pre+ rescaled_math_pre+  rescaled_read_pre +mathxil_z_pre_missing+readxil_z_pre_missing+ oss_dis_pre_zeros+
                incidents_pre_zeros+any_arrests_pre+violent_pre+property_pre+drug_pre, data=d)
    summary(fit)
    #t-test
    output <<- summary(fit)$coefficients["pseudoT",3]
    #pval: #output <<- summary(fit)$coefficients["pseudoT",4]
  }
  runmodel(I)
  return(output)
}

# bootstrapping with 100k replications
ptm <- proc.time()

reps<-100000
results_rescaled_math_post1  <- boot(data=master_dataset_1, statistic=bs_post, strata=master_dataset_1$blocknum,
                R=reps, M="d$rescaled_math_post1", sim="permutation")
results_mathgpa_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$mathgpa_post1", sim="permutation")
results_mathfail_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$mathfail_post1", sim="permutation")

results_rescaled_read_post1 <- boot(data=master_dataset_1, statistic=bs_post, strata=master_dataset_1$blocknum,
                R=reps, M="d$rescaled_read_post1", sim="permutation")
results_nonmathgpa_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$nonmathgpa_post1", sim="permutation")
results_nonmathfail_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$nonmathfail_post1", sim="permutation")

results_incidents_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$incidents_post1", sim="permutation")
results_days_absent_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$days_absent_post1", sim="permutation")
results_oss_dis_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$oss_dis_post1", sim="permutation")

results_violent_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$violent_post1", sim="permutation")
results_property_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$property_post1", sim="permutation")
results_drug_post1 <- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$drug_post1", sim="permutation")
results_other_post1<- boot(data=master_dataset_1, statistic=bs_nopost, strata=master_dataset_1$blocknum,
                R=reps, M="d$other_post1", sim="permutation")

results_complete <- as.data.frame(cbind(results_rescaled_math_post1$t, results_mathgpa_post1$t, results_mathfail_post1$t,
                                        results_rescaled_read_post1$t, results_nonmathgpa_post1$t, results_nonmathfail_post1$t,
                                        results_incidents_post1$t, results_days_absent_post1$t, results_oss_dis_post1$t,
                                        results_violent_post1$t, results_property_post1$t, results_drug_post1$t, results_other_post1$t))
colnames(results_complete) <- c("mathxil_rescaled_post1", "mathgpa_post1", "mathfail_post1",
                                "readxil_rescaled_post1","nonmathgpa_post1", "nonmathfail_post1",
                                "incidents_post1", "absences_post1", "oss_post1",
                                "violent_post1", "property_post1", "drug_post1", "other_post1")
head(results_complete)
proc.time()-ptm

write.csv(results_complete, "/filepath/permutationtests_100k.csv")
