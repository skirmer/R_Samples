# ML Experimentation - Can we predict which CPS grads will complete college? #
# Stephanie Kirmer #
# July 2016 #

# This script uses previously assembled data regarding student academics and predicts their likelihood of completing college (2 year or 4 year diploma)
# Data is proprietary so it can't be shared, sorry!

require(randomForest)
require(caret)
require(dplyr)
require(data.table)
require(ggplot2)
library(pROC)

suppressMessages(source("script that cleaned the data"))

# Data cleaning for modeling #
testdata <- filter(modelfile_schools, !is.na(schlid))

testdata$compss[is.na(testdata$compss)] <- 0
testdata$gpa_act_clg_rdy[is.na(testdata$gpa_act_clg_rdy)] <- 0
testdata$said_ed_flag[is.na(testdata$said_ed_flag)] <- 0
testdata$earliest_4_year[is.na(testdata$earliest_4_year)] <- 0
testdata$earliest_2_year[is.na(testdata$earliest_2_year)] <- 0
testdata$gpa_allmath[is.na(testdata$gpa_allmath)] <- 0
testdata$precalculus[is.na(testdata$precalculus)] <- 0
testdata$calculus[is.na(testdata$calculus)] <- 0
testdata$hs150_200[is.na(testdata$hs150_200)] <- 0
testdata$hs75_150[is.na(testdata$hs75_150)] <- 0
testdata$hs75under[is.na(testdata$hs75under)] <- 0
testdata$hs_200plus[is.na(testdata$hs_200plus)] <- 0
testdata$avg_act[is.na(testdata$avg_act)] <- 0
testdata$avg_gpa[is.na(testdata$avg_gpa)] <- 0
testdata$avgcohortsize[is.na(testdata$avgcohortsize)] <- 0
testdata$earned_anything[is.na(testdata$earned_anything)] <- 0

testdata$race_tree[testdata$white==1] <- 1
testdata$race_tree[testdata$hispanic==1] <- 2
testdata$race_tree[testdata$black==1] <- 3
testdata$race_tree[is.na(testdata$race_tree)] <- 0

testdata$location_tree[testdata$ILonly==1] <- 1
testdata$location_tree[testdata$Chionly==1] <- 2
testdata$location_tree[testdata$NotIL==1] <- 3
testdata$location_tree[is.na(testdata$location_tree)] <- 0

testdata$startplace_tree[testdata$earliest_4_year==1] <- 1
testdata$startplace_tree[testdata$earliest_2_year==1] <- 2
testdata$startplace_tree[is.na(testdata$startplace_tree)] <- 0

testdata$grad_dummy[testdata$GRADUATED == "Y"] <- 1
testdata$grad_dummy[is.na(testdata$GRADUATED)] <- 0

partition <- createDataPartition(y=testdata$earned_anything,
                                 p=.5,
                                 list=F)
training <- testdata[partition,]
testing <- testdata[-partition,]


set.seed(1000)

# ================= Prediction before a student applies to college ================== #####

rf_model_hsonly <- randomForest(factor(grad_dummy) ~ gpa_12th + compss + said_ed_flag + gpa_allmath + max_precalc_or_better +
                                  avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male +
                                  student_match_class_val + race_tree,
                                data = training)

rf_model_hsonly


# The section below that visualizes the importance is borrowed with thanks from Megan L. Ridsal. kaggle.com/mrisdal.
# Get importance
importance    <- importance(rf_model_hsonly)
varImportance <- data.frame(Variables = row.names(importance),
                            Importance = round(importance[ ,'MeanDecreaseGini'],2))

# Create a rank variable based on importance
rankImportance <- varImportance %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))

# Use ggplot2 to visualize the relative importance of variables
ggplot(rankImportance, aes(x = reorder(Variables, Importance),
                           y = Importance, fill = Importance)) +
  geom_bar(stat='identity') +
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'red') +
  labs(x = 'Variables') +
  coord_flip()


# Predict using the test set
prediction <- predict(rf_model_hsonly, testing)

table(prediction)
table(testing$grad_dummy)
table(training$grad_dummy)



# ================= Prediction before a student starts college ================== #####

rf_model_collpick <- randomForest(factor(grad_dummy) ~ gpa_12th + compss + said_ed_flag + startplace_tree + gpa_allmath + max_precalc_or_better +
                                    avgcohortsize + avg_act + avg_gpa + lunch + sped + lep + male + location_tree +
                                    college_class_val+student_match_class_val+ matchplus_new + race_tree,
                                  data = training)


rf_model_collpick


# The section below that visualizes the importance is borrowed with thanks from Megan L. Ridsal. kaggle.com/mrisdal.
# Get importance
importance    <- importance(rf_model_collpick)
varImportance <- data.frame(Variables = row.names(importance),
                            Importance = round(importance[ ,'MeanDecreaseGini'],2))

# Create a rank variable based on importance
rankImportance <- varImportance %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))


# Use ggplot2 to visualize the relative importance of variables
ggplot(rankImportance, aes(x = reorder(Variables, Importance),
                           y = Importance, fill = Importance)) +
  geom_bar(stat='identity') +
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'red') +
  labs(x = 'Variables') +
  coord_flip()


# Predict using the test set
prediction <- predict(rf_model_collpick, testing)


table(prediction)
table(testing$grad_dummy)
table(training$grad_dummy)




# ================= Prediction after a student starts college ================== #####

rf_model_somecoll <- randomForest(factor(grad_dummy) ~ gpa_12th + compss + said_ed_flag + startplace_tree + gpa_allmath + max_precalc_or_better +
                                    avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male + location_tree +
                                    college_class_val+student_match_class_val+ matchplus_new + race_tree + lifetime_college_achievemt,
                                  data = training)

rf_model_somecoll

# The section below that visualizes the importance is borrowed with thanks from Megan L. Ridsal. kaggle.com/mrisdal.
# Get importance
importance    <- importance(rf_model_somecoll)
varImportance <- data.frame(Variables = row.names(importance),
                            Importance = round(importance[ ,'MeanDecreaseGini'],2))

# Create a rank variable based on importance
rankImportance <- varImportance %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))


# Use ggplot2 to visualize the relative importance of variables
ggplot(rankImportance, aes(x = reorder(Variables, Importance),
                           y = Importance, fill = Importance)) +
  geom_bar(stat='identity') +
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'red') +
  labs(x = 'Variables') +
  coord_flip()


# Predict using the test set
prediction <- predict(rf_model_somecoll, testing)

table(prediction)
table(testing$grad_dummy)
table(training$grad_dummy)


# ------------------------------------- ROC - testing the quality of the models ------------------------------------ #
#In plot, specificity = rate of false pos and sensitivity = rate of false neg

gpa_roc <- roc(factor(grad_dummy) ~ gpa_12th,data = training)
act_roc <- roc(factor(grad_dummy) ~ compss,data = training)
ele_roc <- roc(factor(grad_dummy) ~ lifetime_college_achievemt,data = training)


# ==== Model the full kitchen sink approach ==== #####

everything_mod <- glm(grad_dummy~
                        gpa_12th + compss + said_ed_flag + startplace_tree + gpa_allmath + max_precalc_or_better +
                        avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male + location_tree +
                        college_class_val+student_match_class_val+ matchplus_new + race_tree + lifetime_college_achievemt
               , data=training, family=binomial)

#do it with the glm
everything_roc <- roc(everything_mod$y, everything_mod$fitted.values)
plot(everything_roc)


everything_rf <- randomForest(factor(grad_dummy)~
                        gpa_12th + compss + said_ed_flag + startplace_tree + gpa_allmath + max_precalc_or_better +
                        avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male + location_tree +
                        college_class_val+student_match_class_val+ matchplus_new + race_tree + lifetime_college_achievemt
                      , data=training)

# Do it with the random forest
rfModel <- as.vector(predict(everything_rf, testing, type="prob")[,1])
everything_roc <- roc(testing$grad_dummy,rfModel)
plot(everything_roc)

 # ---------------------------------------------------- #####
beforecollege_start_mod <- glm(grad_dummy~
                        gpa_12th + compss + said_ed_flag + startplace_tree + gpa_allmath + max_precalc_or_better +
                        avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male + location_tree +
                        college_class_val+student_match_class_val+ matchplus_new + race_tree
                      , data=training, family=binomial)

#do it with the glm
everything_roc <- roc(beforecollege_start_mod$y, beforecollege_start_mod$fitted.values)
plot(everything_roc)


beforecollege_start_rf <- randomForest(factor(grad_dummy)~
                                         gpa_12th + compss + said_ed_flag + startplace_tree + gpa_allmath + max_precalc_or_better +
                                         avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male + location_tree +
                                         college_class_val+student_match_class_val+ matchplus_new + race_tree
                                       , data=training)

# Do it with the random forest
rfModel <- as.vector(predict(beforecollege_start_rf, testing, type="prob")[,1])
everything_roc <- roc(testing$grad_dummy,rfModel)
plot(everything_roc)

# ---------------------------------------------------- #####
hs_mod <- glm(grad_dummy~
                       gpa_12th + compss + said_ed_flag + gpa_allmath + max_precalc_or_better +
                       avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male +
                        student_match_class_val+ race_tree
                     , data=training, family=binomial)

#do it with the glm
everything_roc <- roc(hs_mod$y, hs_mod$fitted.values)
plot(everything_roc)


hs_rf <- randomForest(factor(grad_dummy)~
                        gpa_12th + compss + said_ed_flag + gpa_allmath + max_precalc_or_better +
                        avgcohortsize+ avg_act + avg_gpa + lunch + sped + lep + male +
                        student_match_class_val+ race_tree
                     , data=training)

# Do it with the random forest
rfModel <- as.vector(predict(hs_rf, testing, type="prob")[,1])
hs_roc <- roc(testing$grad_dummy,rfModel)
plot(hs_roc)



roc.test(everything_roc, hs_roc)

