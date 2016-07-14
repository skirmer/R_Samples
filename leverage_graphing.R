# Script to produce leverage graphs (modeled from Kling 2007) from raw dataset #
# Stephanie Kirmer #
# July 2016 #
# This project is centered around the effect of an intervention in schooling.

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

# CODE STARTS HERE #####
#Data tidying and cleaning has been done before this point - contact me for any questions about that.

itt_for_leverage <- function(outcome, participation1, participation2, grouping, sizefile, weight, number_to_crop, sourcefile, weight_final, omit_type){

  #Assemble formulas for the various models
  itt_step1_form <- formula(paste(outcome, "~", grouping, "+", grouping,":assignment"))
  itt_step2_form <- formula(paste(participation1, "~", grouping, "+", grouping,":assignment"))
  itt_step3optional_form <- formula(paste(participation2, "~", grouping, "+", grouping,":assignment"))
  final_itt_form_1partic <- formula(paste("coeflist_scores~coeflist_partic-1")) #coeflist_partic2
  final_itt_form <- formula(paste("coeflist_scores~coeflist_partic+coeflist_partic2-1")) #coeflist_partic2

  #Get the counts together for weighting later
  blocksizes <<- table(sourcefile[,"blocknum"])
  blocksizes <- as.data.frame(blocksizes)
  schoolsizes <<- table(sourcefile[,"schlid"])
  schoolsizes <- as.data.frame(schoolsizes)

  #adding block or school size for later weighting - block level excludes some blocks but school keeps all
  data_forblocks <- merge(sourcefile, blocksizes, by.x="blocknum", by.y="Var1", all.x=T)
  data_forschools <- merge(sourcefile, schoolsizes, by.x="schlid", by.y="Var1", all.x=T)
  data_forschools$schlid <- as.factor(data_forschools$schlid)
  data_forblocks$schlid <- as.factor(data_forblocks$schlid)

  assign_block <- summarize(group_by(data_forschools, assignment, blocknum),assign_block_count = n())
  assign_school <- summarize(group_by(data_forschools, assignment, schlid),assign_school_count = n())
  school_block <- summarize(group_by(data_forschools, schlid, blocknum),school_block_count = n())
  assign_school_block <- summarize(group_by(data_forschools, assignment, schlid, blocknum),assign_school_block_count = n())
  assign <- summarize(group_by(data_forschools, assignment),assign_count = n())
  allcount <- summarize(group_by(data_forschools),all_count = n())

  data_forschools_withcounts <- merge(data_forschools, assign_block, by=c("assignment", "blocknum"), all.x=T)
  data_forschools_withcounts <- merge(data_forschools_withcounts, assign_school, by=c("assignment", "schlid"), all.x=T)
  data_forschools_withcounts <- merge(data_forschools_withcounts, school_block, by=c("schlid", "blocknum"), all.x=T)
  data_forschools_withcounts <- merge(data_forschools_withcounts, assign_school_block, by=c("assignment", "schlid", "blocknum"), all.x=T)
  data_forschools_withcounts <- merge(data_forschools_withcounts, assign, by=c("assignment"), all.x=T)
  data_forschools_withcounts$all_count <- allcount$all_count

  data_forschools_withcounts<-as.data.table(data_forschools_withcounts, keep.rownames=TRUE)
  data_forschools_withcounts[,mtoweight_calculation:=((assign_count/all_count)/(assign_school_block_count/school_block_count))]
  data_forschools_withcounts[,weight_calculation1:=(assign_count/all_count)]
  data_forschools_withcounts[,weight_calculation2:=(assign_school_block_count/school_block_count)]
  data_forschools_withcounts <- as.data.frame(data_forschools_withcounts)

  if(grouping=="schlid"){
    working_data_file <<- data_forschools_withcounts
    file_string <<- "data_forschools_withcounts"
  }
  if(grouping=="blocknum"){
    working_data_file <<- data_forblocks
    file_string <<- "data_forblocks"
  }

  sample <- summarize(group_by(data_forschools_withcounts, blocknum, assignment, schlid, assign_block_count, assign_school_count, school_block_count, assign_school_block_count, assign_count, all_count, mtoweight_calculation),
                      count = n(),
                      tx_yr1 = sum(treat, na.rm=T),
                      tx_yr2 = sum(treat_y2, na.rm=T),
                      pct_tx_yr1 = tx_yr1/count,
                      pct_tx_yr2 = tx_yr2/count,
                      math_2015_avg = mean(sy15z_mathxil)
  )
  # Layered ITT Work ####

  weight2 <<- eval(parse(text=paste0(file_string, "$", weight)))

  # First ITT - Math Score #
  scoremod <- lm(itt_step1_form, weights=weight2, data=working_data_file)
  coeflist_scores <- coef(scoremod)
  stderr_scores <- coef(summary(scoremod))[,2]
  resid_scores <- resid(scoremod)

  # Second ITT - Year 1 Treatment Participation #
  partic_mod <- lm(itt_step2_form, weights=weight2, data=working_data_file)
  coeflist_partic <- coef(partic_mod)
  stderr_partic <- coef(summary(partic_mod))[,2]
  resid_partic <- resid(partic_mod)

  # Third ITT - Year 2 Treatment Participation #
  partic2_mod <- lm(itt_step3optional_form,weights=weight2, data=working_data_file)
  coeflist_partic2 <- coef(partic2_mod)
  stderr_partic2 <- coef(summary(partic2_mod))[,2]
  resid_partic2 <- resid(partic2_mod)


  # Combine ITT output into a new data frame ####
  scorepartic <- cbind(coeflist_scores, coeflist_partic, coeflist_partic2, stderr_scores, stderr_partic, stderr_partic2)
  scorepartic_resids <- cbind(resid_scores, resid_partic, resid_partic2)

  #Drop the estimates for the individual blocks- we just want the interactions here#
  working_scorepartic <- scorepartic[-c(1:number_to_crop),]
  working_scorepartic <- cbind(sizefile, working_scorepartic)
  working_scorepartic<-as.data.frame(working_scorepartic)

  working_scorepartic<-as.data.table(working_scorepartic, keep.rownames=TRUE)
  working_scorepartic[,stderr_scores_weight:=(1/(working_scorepartic$stderr_scores^2))]
  working_scorepartic[,stderr_partic_weight:=(1/(working_scorepartic$stderr_partic^2))]
  working_scorepartic[,stderr_partic2_weight:=(1/(working_scorepartic$stderr_partic2^2))]
  working_scorepartic <- as.data.frame(working_scorepartic)


  weight_final2 <- eval(parse(text=paste0("working_scorepartic", "$", weight_final)))

  summary(lm(final_itt_form,weights=weight_final2,data=working_scorepartic))
  ittmodel <- lm(final_itt_form,weights=weight_final2,data=working_scorepartic)
  print(ittmodel)

  assign(paste0("ittmodel_", grouping,"_", weight, "_", omit_type), envir=.GlobalEnv, ittmodel) #Create the output with dynamic naming


  leverageform <- formula(paste("J~",  grouping, "+", participation1, "+", participation2, "|", grouping, "+", grouping, ":assignment"))

  suppressMessages(library(formattable))
  suppressMessages(library(knitr))
  suppressMessages(library(data.table))
  suppressMessages(library(dplyr))
  suppressMessages(library(xtable))
  suppressMessages(library(sem))
  suppressMessages(library(AER))
  suppressMessages(library(ivpack))
  fulloutcome <- eval(parse(text=paste0(file_string, "$", outcome)))
  tsls_fn_small <- function(I) {
    J <<- as.matrix(I)
    require(sem)
    model <<- ivreg(leverageform,
                    data=working_data_file)
    totmod <<- summary(model)
    totmod
    robustmod <<- robust.se(model)
    robustmod
  }

  print(tsls_fn_small(fulloutcome))

  assign(paste0("leveragemodel_", grouping,"_", weight, "_", omit_type), envir=.GlobalEnv, model) #Create the output with dynamic naming

}

#run the function
itt_for_leverage(outcome = "sy15z_mathxil", participation1 = "treat_y1_single", participation2 = "treat_y2",
              grouping = "blocknum", sizefile = blocksizes, weight = "precision", number_to_crop = 43,
              sourcefile = math15all, weight_final = "null", "omits")


#Summarize dataset by school/assignment ####
y1bm$schlid_factor <- as.factor(y1bm$schlid)
meansbyschool_3exp <- summarize(group_by(y1bm, schlid_factor, assignment),
                                count = n(),
                                y1_single_mean = mean(treat_y1_single, na.rm=T),
                                y2_single_mean = mean(treat_y2_single, na.rm=T),
                                double_mean = mean(treat_both, na.rm=T),
                                y2double_mean = mean(treat_y2, na.rm=T),

                                y1_daystx_mean = mean(days_treated_y1_zeroes, na.rm=T),
                                y2_daystx_mean = mean(days_treated_y2_zeroes, na.rm=T),
                                both_daystx_mean = mean(days_treated_bothyears_zeroes, na.rm=T),

                                sy13mathmean = mean(sy13z_mathxil, na.rm=T),
                                sy14mathmean = mean(sy14z_mathxil, na.rm=T),
                                sy15mathmean = mean(sy15z_mathxil, na.rm=T)
)

#Summarize dataset by block/assignment ####
meansbyblock_3exp <- summarize(group_by(y1bm, blocknum, assignment),
                               count = n(),
                               y1_single_mean = mean(treat_y1_single, na.rm=T),
                               y2_single_mean = mean(treat_y2_single, na.rm=T),
                               double_mean = mean(treat_both, na.rm=T),
                               y2double_mean = mean(treat_y2, na.rm=T),

                               y1_daystx_mean = mean(days_treated_y1_zeroes, na.rm=T),
                               y2_daystx_mean = mean(days_treated_y2_zeroes, na.rm=T),
                               both_daystx_mean = mean(days_treated_bothyears_zeroes, na.rm=T),

                               sy13mathmean = mean(sy13z_mathxil, na.rm=T),
                               sy14mathmean = mean(sy14z_mathxil, na.rm=T),
                               sy15mathmean = mean(sy15z_mathxil, na.rm=T)
)

#Summarize dataset by school/NOT assignment ####
meansbyschool_3_noassign <- summarize(group_by(y1bm, schlid_factor),
                                      y1_single_mean_sch = mean(treat_y1_single, na.rm=T),
                                      y2_single_mean_sch = mean(treat_y2_single, na.rm=T),
                                      double_mean_sch = mean(treat_both, na.rm=T),
                                      y2double_mean_sch = mean(treat_y2, na.rm=T),

                                      y1_daystx_mean_sch = mean(days_treated_y1_zeroes, na.rm=T),
                                      y2_daystx_mean_sch = mean(days_treated_y2_zeroes, na.rm=T),
                                      both_daystx_mean_sch = mean(days_treated_bothyears_zeroes, na.rm=T),

                                      sy13mathmean_sch = mean(sy13z_mathxil, na.rm=T),
                                      sy14mathmean_sch = mean(sy14z_mathxil, na.rm=T),
                                      sy15mathmean_sch = mean(sy15z_mathxil, na.rm=T)
)

#Summarize dataset by block/NOT assignment ####

meansbyblock_3_noassign <- summarize(group_by(y1bm, blocknum),
                                     y1_single_mean_bk = mean(treat_y1_single, na.rm=T),
                                     y2_single_mean_bk = mean(treat_y2_single, na.rm=T),
                                     double_mean_bk = mean(treat_both, na.rm=T),

                                     y1_daystx_mean_bk = mean(days_treated_y1_zeroes, na.rm=T),
                                     y2_daystx_mean_bk = mean(days_treated_y2_zeroes, na.rm=T),
                                     both_daystx_mean_bk = mean(days_treated_bothyears_zeroes, na.rm=T),

                                     sy13mathmean_bk = mean(sy13z_mathxil, na.rm=T),
                                     sy14mathmean_bk = mean(sy14z_mathxil, na.rm=T),
                                     sy15mathmean_bk = mean(sy15z_mathxil, na.rm=T)
)

# Combine the summarized files to allow for decentering means ####
merged_means <- merge(meansbyblock_3exp, meansbyblock_3_noassign, by="blocknum", all.x=T)
merged_means_sch <- merge(meansbyschool_3exp, meansbyschool_3_noassign, by="schlid_factor", all.x=T)

merged_means$decentered_2014_math <- merged_means$sy14mathmean - merged_means$sy14mathmean_bk
merged_means$decentered_2015_math <- merged_means$sy15mathmean - merged_means$sy15mathmean_bk

merged_means_sch$decentered_2014_math <- merged_means_sch$sy14mathmean - merged_means_sch$sy14mathmean_sch
merged_means_sch$decentered_2015_math <- merged_means_sch$sy15mathmean - merged_means_sch$sy15mathmean_sch


# Block Level leverage Leverages - precision, omit ####
model <- coef(leveragemodel_blocknum_precision_omits)

block_all_prcs_y1 <- ggplot(merged_means, aes(y1_single_mean,decentered_2015_math, label=blocknum))+
  theme(panel.background=element_rect(fill="white", color="black"), # Set up the background
        text=element_text(family="Times"), # Default font is serif
        plot.title=element_text(size=rel(1.5), hjust=0), # Make the title left justified and a bit bigger than other text
        axis.text.y=element_text(size=14), # Y axis text size
        axis.text.x=element_text(size=14), # X axis text size
        axis.title.y=element_text(size=16), # Y axis title text size
        axis.title.x=element_text(size=16), # X axis title text size
        legend.title=element_text(size=14), # Legend title text size
        legend.text=element_text(size=13), # Legend text size
        legend.justification=c(1,0), #Legend location pt 1
        legend.position=c(1,0), # Legend location pt 2
        legend.background=element_rect(color="black"))+ # Outline the legend
  geom_point(aes(fill=as.factor(assignment), size=count), stroke=.5, color="black", pch=21)+ #Points are filled according to assignment, sized according to size, outlined in black
  guides(fill=guide_legend(override.aes=list(size=10)))+ # Make the legend color points a big bigger for viewing
  geom_text_repel(aes(label=blocknum),nudge_x=0.015, nudge_y=-0.01, size=4.5, family="Times")+ # Scoot the labels around inside the plot, with lead lines
  scale_fill_manual(values=c("#f0f0f0","#969696"), name="Assignment", labels=c("Control", "Treatment"))+ # What the colors should be and titles in the legend
  scale_size(range=c(0,20), "# of Tx Students/Block/Assignment", guide="none")+ # Range of point sizes and the naming if there were a legend, but guide=none means there's no legend
  #geom_text(check_overlap=FALSE, nudge_x=0.02, nudge_y=0.005)+ #Hidden to allow geom_text_repel to work
  geom_abline(intercept=0, slope=model["treat_y1_single"])+ # Add the fitted line using the slope from our model
  labs(title="Figure X: Y1 Participation Mean by Year 2 Math Score (2015) by Block", y="Y2 Math", x="Y1 Single Dose Mean") # Label the actual plot
block_all_prcs_y1


model <- coef(leveragemodel_blocknum_precision_noomits)

block_all_prcs_y2 <- ggplot(merged_means, aes(y2double_mean,decentered_2015_math, label=blocknum))+
  theme(panel.background=element_rect(fill="white", color="black"), # Set up the background
        text=element_text(family="Times"), # Default font is serif
        plot.title=element_text(size=rel(1.5), hjust=0), # Make the title left justified and a bit bigger than other text
        axis.text.y=element_text(size=14), # Y axis text size
        axis.text.x=element_text(size=14), # X axis text size
        axis.title.y=element_text(size=16), # Y axis title text size
        axis.title.x=element_text(size=16), # X axis title text size
        legend.title=element_text(size=14), # Legend title text size
        legend.text=element_text(size=13), # Legend text size
        legend.justification=c(1,0), #Legend location pt 1
        legend.position=c(1,0), # Legend location pt 2
        legend.background=element_rect(color="black"))+ # Outline the legend
  geom_point(aes(fill=as.factor(assignment), size=count), stroke=.5, color="black", pch=21)+ #Points are filled according to assignment, sized according to size, outlined in black
  guides(fill=guide_legend(override.aes=list(size=10)))+ # Make the legend color points a big bigger for viewing
  geom_text_repel(aes(label=blocknum),nudge_x=0.015, nudge_y=-0.01, size=4.5, family="Times")+ # Scoot the labels around inside the plot, with lead lines
  scale_fill_manual(values=c("#f0f0f0","#969696"), name="Assignment", labels=c("Control", "Treatment"))+ # What the colors should be and titles in the legend
  scale_size(range=c(0,20), "# of Tx Students/Block/Assignment", guide="none")+ # Range of point sizes and the naming if there were a legend, but guide=none means there's no legend
  #geom_text(check_overlap=FALSE, nudge_x=0.02, nudge_y=0.005)+ #Hidden to allow geom_text_repel to work
  geom_abline(intercept=0, slope=model["treat_y2"])+ # Add the fitted line using the slope from our model
  labs(title="Figure X: Y2 Participation Mean by Year 2 Math Score (2015) by Block", y="Y2 Math", x="Y2 Any Dose Mean") # Label the actual plot

block_all_prcs_y2

