#####################################################################
# Response Monitoring Theta-Band Activities across Emotional Contexts in Schizophrenia- and Bipolar-Spectrum Disorders
# Suzuki, Menkes, et al.
#
# Written by Takakuni Suzuki and adapted by Margo W. Menkes
#####################################################################

rm(list=ls())

library(psych)
library(plyr)
library(data.table)

#Stim markers: 1, 2 
#     1 = congruent 
#     2 = incongruent 
#     
#   Response markers: 3, 4, 9, 10
#     3 or 4 = correct 
#     9 or 10 = incorrect

############################ 
############################ FARR
############################
setwd('')


part_list_FARR <- read.table("./Output_Auto_Combined/EEG_FARR_Data_Resp_Combined_2024-07-06.csv", header = TRUE, sep = ",")
setwd("./Output_Auto_Combined/FARR/preprocess_2_behavioral_dbl") 

file_list_FARR <- list.files()

Beh_FARR <- do.call("rbind", lapply(file_list_FARR, function(x) {
  dat <- read.table(x, header=TRUE, sep="\t", fill = TRUE)
  dat$ID <- tools::file_path_sans_ext(basename(x))
  dat
}))

Beh_FARR$ID <- substr(Beh_FARR$ID, 1, 13) #to keep only TNA_FARR_#### in file name
Beh1_FARR <- subset(Beh_FARR, Beh_FARR$ID %in% part_list_FARR$ID)
#Remove irrelevant trials
Beh1_FARR <- Beh1_FARR[!Beh1_FARR$type == "boundary", ]
Beh1_FARR <- Beh1_FARR[!Beh1_FARR$type == "S 62", ]
Beh1_FARR <- Beh1_FARR[(Beh1_FARR$type == "S  1" | Beh1_FARR$type == "S  2" | Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4"| Beh1_FARR$type == "S  9"|Beh1_FARR$type == "S 10"), ]

#Check number of observations for each paticipant
Beh1_FARR$EC <- NA
Beh1_FARR$EC=ifelse((Beh1_FARR$type=="S  9"), 0, 
               ifelse((Beh1_FARR$type=="S 10"), 0, 
                      ifelse((Beh1_FARR$type=="S  3"), 1, 
                             ifelse((Beh1_FARR$type=="S  4"), 1, NA))))
# Calculate RT
Beh1_FARR$latency0 <- shift(Beh1_FARR$latency, n = 1L)

Beh1_FARR$RT <- NA
Beh1_FARR$RT = ifelse((Beh1_FARR$type=="S  9"), Beh1_FARR$latency - Beh1_FARR$latency0, 
                 ifelse((Beh1_FARR$type=="S 10"), Beh1_FARR$latency - Beh1_FARR$latency0, 
                        ifelse((Beh1_FARR$type=="S  3"), Beh1_FARR$latency - Beh1_FARR$latency0, 
                               ifelse((Beh1_FARR$type=="S  4"), Beh1_FARR$latency - Beh1_FARR$latency0, NA))))

Beh1_FARR$prev <- shift(Beh1_FARR$type, n = 1L)
Beh1_FARR$prev_n <- shift(Beh1_FARR$ID, n = 1L)
Beh1_FARR$RT <- ifelse(Beh1_FARR$ID == Beh1_FARR$prev_n, Beh1_FARR$RT, NA) #calculating RT using all trials


Beh1_FARR_MAD <- aggregate(RT~ID, Beh1_FARR, function(x) c(median = median(x), mad = mad(x))) 
Beh1_FARR_MAD$median_min <- Beh1_FARR_MAD$RT[,"median"] - 3*Beh1_FARR_MAD$RT[,"mad"]
Beh1_FARR_MAD$median_max <- Beh1_FARR_MAD$RT[,"median"] + 3*Beh1_FARR_MAD$RT[,"mad"]
Beh1_FARR <- merge(Beh1_FARR, Beh1_FARR_MAD[,c("ID", "median_min", "median_max")], all.x = TRUE, by = "ID")

Beh1_FARR$RT_trimmed <- ifelse(Beh1_FARR$RT > .1 & Beh1_FARR$RT < 1.5 & Beh1_FARR$RT > Beh1_FARR$median_min & Beh1_FARR$RT < Beh1_FARR$median_max, Beh1_FARR$RT, NA) # calc RT using only trials w/ BOTH an RT btw 100 - 1500 ms & w/in median min & max for this person 


# Various RT based calculations 
Beh1_FARR$PES1 <- shift(Beh1_FARR$type, n = 2L) # shift event marker down 
Beh1_FARR$PES1_n <- shift(Beh1_FARR$ID, n = 2L) # also shift IDs for the row

Beh1_FARR$PESneg1 <- shift(Beh1_FARR$type, n = -2L) # shift event marker up (indicate what is response for after this trial)
Beh1_FARR$PESneg1_n <- shift(Beh1_FARR$ID, n = -2L) # also shift IDs 

#Remove any overlap from previous participant
Beh1_FARR$PES1 <- ifelse(Beh1_FARR$ID == Beh1_FARR$PES1_n, Beh1_FARR$PES1, NA)
Beh1_FARR$PESneg1 <- ifelse(Beh1_FARR$ID == Beh1_FARR$PESneg1_n, Beh1_FARR$PESneg1, NA)

#Keep correct RTs
Beh1_FARR$CorrectRT <- ifelse((Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4"), Beh1_FARR$RT_trimmed, NA)

#Keep Error RTs
Beh1_FARR$ErrorRT <- ifelse((Beh1_FARR$type == "S  9" | Beh1_FARR$type == "S 10"), Beh1_FARR$RT_trimmed, NA)

#Keep correct after correct number RTs
Beh1_FARR$CAfterCRT <- ifelse(((Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4")), Beh1_FARR$RT_trimmed, NA)

#Keep error after correct number RTs
Beh1_FARR$EAfterCRT <- ifelse(((Beh1_FARR$type == "S  9" | Beh1_FARR$type == "S 10") & (Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4")), Beh1_FARR$RT_trimmed, NA)

#Find post-error correct RTs (this trial has to be correct)
Beh1_FARR$CAfterERT <- ifelse(((Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$PES1 == "S  9" | Beh1_FARR$PES1 == "S 10")), Beh1_FARR$RT_trimmed, NA)

#post Error correct trials
Beh1_FARR$PostError_inc <- ifelse(((Beh1_FARR$PES1 == "S  9" | Beh1_FARR$PES1 == "S 10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  2")), Beh1_FARR$RT_trimmed, NA)
Beh1_FARR$PostError_con <- ifelse(((Beh1_FARR$PES1 == "S  9" | Beh1_FARR$PES1 == "S  10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  1")), Beh1_FARR$RT_trimmed, NA)

#Pre Error correct trials
Beh1_FARR$Pre_Error_inc <- ifelse(((Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4") & (Beh1_FARR$PESneg1 == "S  9" | Beh1_FARR$PESneg1 == "S 10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  2")), Beh1_FARR$RT_trimmed, NA)
Beh1_FARR$Pre_Error_con <- ifelse(((Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4") & (Beh1_FARR$PESneg1 == "S  9" | Beh1_FARR$PESneg1 == "S 10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  1")), Beh1_FARR$RT_trimmed, NA)

#post Error accuracy 
Beh1_FARR$PostError_inc_EC <- ifelse(((Beh1_FARR$PES1 == "S  9" | Beh1_FARR$PES1 == "S 10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  2")), 1, 
                                ifelse(((Beh1_FARR$PES1 == "S  9" | Beh1_FARR$PES1 == "S 10") & (Beh1_FARR$type == "S  9" | Beh1_FARR$type == "S 10") & (Beh1_FARR$prev == "S  2")), 0, NA))
Beh1_FARR$PostError_con_EC <- ifelse(((Beh1_FARR$PES1 == "S  9" | Beh1_FARR$PES1 == "S 10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  1")), 1,
                                ifelse(((Beh1_FARR$PES1 == "S  9" | Beh1_FARR$PES1 == "S 10") & (Beh1_FARR$type == "S  9" | Beh1_FARR$type == "S 10") & (Beh1_FARR$prev == "S  1")), 0, NA))

#Pre Error accuracy
Beh1_FARR$Pre_Error_inc_EC <- ifelse(((Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4") & (Beh1_FARR$PESneg1 == "S  9" | Beh1_FARR$PESneg1 == "S 10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  2")), 1,
                             ifelse(((Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4") & (Beh1_FARR$PESneg1 == "S  9" | Beh1_FARR$PESneg1 == "S 10") & (Beh1_FARR$type == "S  9" | Beh1_FARR$type == "S 10") & (Beh1_FARR$prev == "S  2")), 0, NA))
Beh1_FARR$Pre_Error_con_EC <- ifelse(((Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4") & (Beh1_FARR$PESneg1 == "S  9" | Beh1_FARR$PESneg1 == "S 10") & (Beh1_FARR$type == "S  3" | Beh1_FARR$type == "S  4") & (Beh1_FARR$prev == "S  1")), 1,
                             ifelse(((Beh1_FARR$PES1 == "S  3" | Beh1_FARR$PES1 == "S  4") & (Beh1_FARR$PESneg1 == "S  9" | Beh1_FARR$PESneg1 == "S 10") & (Beh1_FARR$type == "S  9" | Beh1_FARR$type == "S 10") & (Beh1_FARR$prev == "S  1")), 0, NA))

#Congruent number RTs
Beh1_FARR$CongruentRT <- ifelse((Beh1_FARR$prev == "S  1"), Beh1_FARR$RT_trimmed, NA)

#Incongruent number RTs

Beh1_FARR$IncongruentRT <- ifelse((Beh1_FARR$prev == "S  2"), Beh1_FARR$RT_trimmed, NA)

#Coding for Congruent vs. Incongruent and Correct vs. Error
Beh1_FARR$Cong_Correct <- ifelse((Beh1_FARR$prev == "S  1"), ifelse(Beh1_FARR$EC == 1, 1, 0), NA)
Beh1_FARR$Incong_Correct <- ifelse((Beh1_FARR$prev == "S  2"), ifelse(Beh1_FARR$EC == 1, 1, 0), NA)

#Calculate averages
BehAve_FARR <- aggregate(Beh1_FARR[, c("EC", "RT", "RT_trimmed", "CorrectRT", "ErrorRT", "CAfterCRT", "EAfterCRT", "CAfterERT", 
                                       "PostError_inc", "PostError_con", "Pre_Error_inc", "Pre_Error_con",
                                       "CongruentRT", "IncongruentRT", "Cong_Correct", "Incong_Correct")], list(Beh1_FARR$ID), mean, na.rm=TRUE)
BehAve_FARR$PES_Trad <- BehAve_FARR$CAfterERT - BehAve_FARR$CAfterCRT
BehAve_FARR$PES_Robust <- (BehAve_FARR$PostError_inc + BehAve_FARR$PostError_con)/2 - (BehAve_FARR$Pre_Error_inc + BehAve_FARR$Pre_Error_con)/2


############################ 
############################ FNEG 
############################

part_list_FNEG <- read.table("./Output_Auto_Combined/EEG_FNEG_Data_Resp_Combined_2024-07-06.csv", header = TRUE, sep = ",")
setwd("./Output_Auto_Combined/FNEG/preprocess_2_behavioral_dbl") 

file_list_FNEG <- list.files()

Beh_FNEG <- do.call("rbind", lapply(file_list_FNEG, function(x) {
  dat <- read.table(x, header=TRUE, sep="\t", fill = TRUE)
  dat$ID <- tools::file_path_sans_ext(basename(x))
  dat
}))

Beh_FNEG$ID <- substr(Beh_FNEG$ID, 1, 13) #to keep only TNA_FNEG_#### in file name
Beh1_FNEG <- subset(Beh_FNEG, Beh_FNEG$ID %in% part_list_FNEG$ID)
#Remove irrelevant trials
Beh1_FNEG <- Beh1_FNEG[!Beh1_FNEG$type == "boundary", ]
Beh1_FNEG <- Beh1_FNEG[!Beh1_FNEG$type == "S 62", ]
Beh1_FNEG <- Beh1_FNEG[(Beh1_FNEG$type == "S  1" | Beh1_FNEG$type == "S  2" | Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4"| Beh1_FNEG$type == "S  9"|Beh1_FNEG$type == "S 10"), ]

#Check number of observations for each paticipant

#Make column for correct=1 and error=0
Beh1_FNEG$EC <- NA
Beh1_FNEG$EC=ifelse((Beh1_FNEG$type=="S  9"), 0, 
                    ifelse((Beh1_FNEG$type=="S 10"), 0, 
                           ifelse((Beh1_FNEG$type=="S  3"), 1, 
                                  ifelse((Beh1_FNEG$type=="S  4"), 1, NA))))
# Calculate RT
library(data.table)
Beh1_FNEG$latency0 <- shift(Beh1_FNEG$latency, n = 1L)

Beh1_FNEG$RT <- NA
Beh1_FNEG$RT = ifelse((Beh1_FNEG$type=="S  9"), Beh1_FNEG$latency - Beh1_FNEG$latency0, 
                      ifelse((Beh1_FNEG$type=="S 10"), Beh1_FNEG$latency - Beh1_FNEG$latency0, 
                             ifelse((Beh1_FNEG$type=="S  3"), Beh1_FNEG$latency - Beh1_FNEG$latency0, 
                                    ifelse((Beh1_FNEG$type=="S  4"), Beh1_FNEG$latency - Beh1_FNEG$latency0, NA))))

Beh1_FNEG$prev <- shift(Beh1_FNEG$type, n = 1L)
Beh1_FNEG$prev_n <- shift(Beh1_FNEG$ID, n = 1L)
Beh1_FNEG$RT <- ifelse(Beh1_FNEG$ID == Beh1_FNEG$prev_n, Beh1_FNEG$RT, NA) #calculating RT using all trials

Beh1_FNEG_MAD <- aggregate(RT~ID, Beh1_FNEG, function(x) c(median = median(x), mad = mad(x)))
Beh1_FNEG_MAD$median_min <- Beh1_FNEG_MAD$RT[,"median"] - 3*Beh1_FNEG_MAD$RT[,"mad"]
Beh1_FNEG_MAD$median_max <- Beh1_FNEG_MAD$RT[,"median"] + 3*Beh1_FNEG_MAD$RT[,"mad"]
Beh1_FNEG <- merge(Beh1_FNEG, Beh1_FNEG_MAD[,c("ID", "median_min", "median_max")], all.x = TRUE, by = "ID")

Beh1_FNEG$RT_trimmed <- ifelse(Beh1_FNEG$RT > .1 & Beh1_FNEG$RT < 1.5 & Beh1_FNEG$RT > Beh1_FNEG$median_min & Beh1_FNEG$RT < Beh1_FNEG$median_max, Beh1_FNEG$RT, NA) # calc RT using only trials w/ BOTH an RT btw 100 - 1500 ms & w/in median min & max for this person 


# Various RT based calculations 
Beh1_FNEG$PES1 <- shift(Beh1_FNEG$type, n = 2L)
Beh1_FNEG$PES1_n <- shift(Beh1_FNEG$ID, n = 2L)

Beh1_FNEG$PESneg1 <- shift(Beh1_FNEG$type, n = -2L)
Beh1_FNEG$PESneg1_n <- shift(Beh1_FNEG$ID, n = -2L)

#Remove any overlap from previous participant
Beh1_FNEG$PES1 <- ifelse(Beh1_FNEG$ID == Beh1_FNEG$PES1_n, Beh1_FNEG$PES1, NA)
Beh1_FNEG$PESneg1 <- ifelse(Beh1_FNEG$ID == Beh1_FNEG$PESneg1_n, Beh1_FNEG$PESneg1, NA)

#Keep correct RTs
Beh1_FNEG$CorrectRT <- ifelse((Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4"), Beh1_FNEG$RT_trimmed, NA)

#Keep Error RTs
Beh1_FNEG$ErrorRT <- ifelse((Beh1_FNEG$type == "S  9" | Beh1_FNEG$type == "S 10"), Beh1_FNEG$RT_trimmed, NA)

#Keep correct after correct number RTs
Beh1_FNEG$CAfterCRT <- ifelse(((Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4")), Beh1_FNEG$RT_trimmed, NA)

#Keep error after correct number RTs
Beh1_FNEG$EAfterCRT <- ifelse(((Beh1_FNEG$type == "S  9" | Beh1_FNEG$type == "S 10") & (Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4")), Beh1_FNEG$RT_trimmed, NA)

#Find post-error correct RTs (this trial has to be correct)
Beh1_FNEG$CAfterERT <- ifelse(((Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$PES1 == "S  9" | Beh1_FNEG$PES1 == "S 10")), Beh1_FNEG$RT_trimmed, NA)

#post Error correct trials
Beh1_FNEG$PostError_inc <- ifelse(((Beh1_FNEG$PES1 == "S  9" | Beh1_FNEG$PES1 == "S 10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  2")), Beh1_FNEG$RT_trimmed, NA)
Beh1_FNEG$PostError_con <- ifelse(((Beh1_FNEG$PES1 == "S  9" | Beh1_FNEG$PES1 == "S  10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  1")), Beh1_FNEG$RT_trimmed, NA)

#Pre Error correct trials
Beh1_FNEG$Pre_Error_inc <- ifelse(((Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4") & (Beh1_FNEG$PESneg1 == "S  9" | Beh1_FNEG$PESneg1 == "S 10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  2")), Beh1_FNEG$RT_trimmed, NA)
Beh1_FNEG$Pre_Error_con <- ifelse(((Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4") & (Beh1_FNEG$PESneg1 == "S  9" | Beh1_FNEG$PESneg1 == "S 10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  1")), Beh1_FNEG$RT_trimmed, NA)

#post Error accuracy 
Beh1_FNEG$PostError_inc_EC <- ifelse(((Beh1_FNEG$PES1 == "S  9" | Beh1_FNEG$PES1 == "S 10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  2")), 1, 
                                     ifelse(((Beh1_FNEG$PES1 == "S  9" | Beh1_FNEG$PES1 == "S 10") & (Beh1_FNEG$type == "S  9" | Beh1_FNEG$type == "S 10") & (Beh1_FNEG$prev == "S  2")), 0, NA))
Beh1_FNEG$PostError_con_EC <- ifelse(((Beh1_FNEG$PES1 == "S  9" | Beh1_FNEG$PES1 == "S 10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  1")), 1,
                                     ifelse(((Beh1_FNEG$PES1 == "S  9" | Beh1_FNEG$PES1 == "S 10") & (Beh1_FNEG$type == "S  9" | Beh1_FNEG$type == "S 10") & (Beh1_FNEG$prev == "S  1")), 0, NA))

#Pre Error accuracy
Beh1_FNEG$Pre_Error_inc_EC <- ifelse(((Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4") & (Beh1_FNEG$PESneg1 == "S  9" | Beh1_FNEG$PESneg1 == "S 10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  2")), 1,
                                     ifelse(((Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4") & (Beh1_FNEG$PESneg1 == "S  9" | Beh1_FNEG$PESneg1 == "S 10") & (Beh1_FNEG$type == "S  9" | Beh1_FNEG$type == "S 10") & (Beh1_FNEG$prev == "S  2")), 0, NA))
Beh1_FNEG$Pre_Error_con_EC <- ifelse(((Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4") & (Beh1_FNEG$PESneg1 == "S  9" | Beh1_FNEG$PESneg1 == "S 10") & (Beh1_FNEG$type == "S  3" | Beh1_FNEG$type == "S  4") & (Beh1_FNEG$prev == "S  1")), 1,
                                     ifelse(((Beh1_FNEG$PES1 == "S  3" | Beh1_FNEG$PES1 == "S  4") & (Beh1_FNEG$PESneg1 == "S  9" | Beh1_FNEG$PESneg1 == "S 10") & (Beh1_FNEG$type == "S  9" | Beh1_FNEG$type == "S 10") & (Beh1_FNEG$prev == "S  1")), 0, NA))

#Congruent number RTs
Beh1_FNEG$CongruentRT <- ifelse((Beh1_FNEG$prev == "S  1"), Beh1_FNEG$RT_trimmed, NA)

#Incongruent number RTs
Beh1_FNEG$IncongruentRT <- ifelse((Beh1_FNEG$prev == "S  2"), Beh1_FNEG$RT_trimmed, NA)

#Coding for Congruent vs. Incongruent and Correct vs. Error
Beh1_FNEG$Cong_Correct <- ifelse((Beh1_FNEG$prev == "S  1"), ifelse(Beh1_FNEG$EC == 1, 1, 0), NA)
Beh1_FNEG$Incong_Correct <- ifelse((Beh1_FNEG$prev == "S  2"), ifelse(Beh1_FNEG$EC == 1, 1, 0), NA)


#Calculate averages
BehAve_FNEG <- aggregate(Beh1_FNEG[, c("EC", "RT", "RT_trimmed", "CorrectRT", "ErrorRT", "CAfterCRT", "EAfterCRT", "CAfterERT", 
                                       "PostError_inc", "PostError_con", "Pre_Error_inc", "Pre_Error_con",
                                       "CongruentRT", "IncongruentRT", "Cong_Correct", "Incong_Correct")], list(Beh1_FNEG$ID), mean, na.rm=TRUE)
BehAve_FNEG$PES_Trad <- BehAve_FNEG$CAfterERT - BehAve_FNEG$CAfterCRT
BehAve_FNEG$PES_Robust <- (BehAve_FNEG$PostError_inc + BehAve_FNEG$PostError_con)/2 - (BehAve_FNEG$Pre_Error_inc + BehAve_FNEG$Pre_Error_con)/2


############################ 
############################ FPOS
############################

part_list_FPOS <- read.table("./Output_Auto_Combined/EEG_FPOS_Data_Resp_Combined_2024-07-06.csv", header = TRUE, sep = ",")
setwd("./Output_Auto_Combined/FPOS/preprocess_2_behavioral_dbl") 

file_list_FPOS <- list.files()

Beh_FPOS <- do.call("rbind", lapply(file_list_FPOS, function(x) {
  dat <- read.table(x, header=TRUE, sep="\t", fill = TRUE)
  dat$ID <- tools::file_path_sans_ext(basename(x))
  dat
}))

Beh_FPOS$ID <- substr(Beh_FPOS$ID, 1, 13) #to keep only TNA_FNEG_#### in file name
Beh1_FPOS <- subset(Beh_FPOS, Beh_FPOS$ID %in% part_list_FPOS$ID)
#Remove irrelevant trials
Beh1_FPOS <- Beh1_FPOS[!Beh1_FPOS$type == "boundary", ]
Beh1_FPOS <- Beh1_FPOS[!Beh1_FPOS$type == "S 62", ]
Beh1_FPOS <- Beh1_FPOS[(Beh1_FPOS$type == "S  1" | Beh1_FPOS$type == "S  2" | Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4"| Beh1_FPOS$type == "S  9"|Beh1_FPOS$type == "S 10"), ]

#Check number of observations for each paticipant

#Make column for correct=1 and error=0
Beh1_FPOS$EC <- NA
Beh1_FPOS$EC=ifelse((Beh1_FPOS$type=="S  9"), 0, 
                    ifelse((Beh1_FPOS$type=="S 10"), 0, 
                           ifelse((Beh1_FPOS$type=="S  3"), 1, 
                                  ifelse((Beh1_FPOS$type=="S  4"), 1, NA))))
# Calculate RT
Beh1_FPOS$latency0 <- shift(Beh1_FPOS$latency, n = 1L)

Beh1_FPOS$RT <- NA
Beh1_FPOS$RT = ifelse((Beh1_FPOS$type=="S  9"), Beh1_FPOS$latency - Beh1_FPOS$latency0, 
                      ifelse((Beh1_FPOS$type=="S 10"), Beh1_FPOS$latency - Beh1_FPOS$latency0, 
                             ifelse((Beh1_FPOS$type=="S  3"), Beh1_FPOS$latency - Beh1_FPOS$latency0, 
                                    ifelse((Beh1_FPOS$type=="S  4"), Beh1_FPOS$latency - Beh1_FPOS$latency0, NA))))

Beh1_FPOS$prev <- shift(Beh1_FPOS$type, n = 1L)
Beh1_FPOS$prev_n <- shift(Beh1_FPOS$ID, n = 1L)
Beh1_FPOS$RT <- ifelse(Beh1_FPOS$ID == Beh1_FPOS$prev_n, Beh1_FPOS$RT, NA) #calculating RT using all trials


Beh1_FPOS_MAD <- aggregate(RT~ID, Beh1_FPOS, function(x) c(median = median(x), mad = mad(x)))
Beh1_FPOS_MAD$median_min <- Beh1_FPOS_MAD$RT[,"median"] - 3*Beh1_FPOS_MAD$RT[,"mad"]
Beh1_FPOS_MAD$median_max <- Beh1_FPOS_MAD$RT[,"median"] + 3*Beh1_FPOS_MAD$RT[,"mad"]
Beh1_FPOS <- merge(Beh1_FPOS, Beh1_FPOS_MAD[,c("ID", "median_min", "median_max")], all.x = TRUE, by = "ID")

Beh1_FPOS$RT_trimmed <- ifelse(Beh1_FPOS$RT > .1 & Beh1_FPOS$RT < 1.5 & Beh1_FPOS$RT > Beh1_FPOS$median_min & Beh1_FPOS$RT < Beh1_FPOS$median_max, Beh1_FPOS$RT, NA) # calc RT using only trials w/ BOTH an RT btw 100 - 1500 ms & w/in median min & max for this person 

# Various RT based calculations 
Beh1_FPOS$PES1 <- shift(Beh1_FPOS$type, n = 2L)
Beh1_FPOS$PES1_n <- shift(Beh1_FPOS$ID, n = 2L)

Beh1_FPOS$PESneg1 <- shift(Beh1_FPOS$type, n = -2L)
Beh1_FPOS$PESneg1_n <- shift(Beh1_FPOS$ID, n = -2L)

#Remove any overlap from previous participant
Beh1_FPOS$PES1 <- ifelse(Beh1_FPOS$ID == Beh1_FPOS$PES1_n, Beh1_FPOS$PES1, NA)
Beh1_FPOS$PESneg1 <- ifelse(Beh1_FPOS$ID == Beh1_FPOS$PESneg1_n, Beh1_FPOS$PESneg1, NA)

#Keep correct RTs
Beh1_FPOS$CorrectRT <- ifelse((Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4"), Beh1_FPOS$RT_trimmed, NA)

#Keep Error RTs
Beh1_FPOS$ErrorRT <- ifelse((Beh1_FPOS$type == "S  9" | Beh1_FPOS$type == "S 10"), Beh1_FPOS$RT_trimmed, NA)

#Keep correct after correct number RTs
Beh1_FPOS$CAfterCRT <- ifelse(((Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4")), Beh1_FPOS$RT_trimmed, NA)

#Keep error after correct number RTs
Beh1_FPOS$EAfterCRT <- ifelse(((Beh1_FPOS$type == "S  9" | Beh1_FPOS$type == "S 10") & (Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4")), Beh1_FPOS$RT_trimmed, NA)

#Find post-error correct RTs (this trial has to be correct)
Beh1_FPOS$CAfterERT <- ifelse(((Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$PES1 == "S  9" | Beh1_FPOS$PES1 == "S 10")), Beh1_FPOS$RT_trimmed, NA)

#post Error correct trials
Beh1_FPOS$PostError_inc <- ifelse(((Beh1_FPOS$PES1 == "S  9" | Beh1_FPOS$PES1 == "S 10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  2")), Beh1_FPOS$RT_trimmed, NA)
Beh1_FPOS$PostError_con <- ifelse(((Beh1_FPOS$PES1 == "S  9" | Beh1_FPOS$PES1 == "S  10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  1")), Beh1_FPOS$RT_trimmed, NA)

#Pre Error correct trials
Beh1_FPOS$Pre_Error_inc <- ifelse(((Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4") & (Beh1_FPOS$PESneg1 == "S  9" | Beh1_FPOS$PESneg1 == "S 10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  2")), Beh1_FPOS$RT_trimmed, NA)
Beh1_FPOS$Pre_Error_con <- ifelse(((Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4") & (Beh1_FPOS$PESneg1 == "S  9" | Beh1_FPOS$PESneg1 == "S 10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  1")), Beh1_FPOS$RT_trimmed, NA)

#post Error accuracy 
Beh1_FPOS$PostError_inc_EC <- ifelse(((Beh1_FPOS$PES1 == "S  9" | Beh1_FPOS$PES1 == "S 10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  2")), 1, 
                                     ifelse(((Beh1_FPOS$PES1 == "S  9" | Beh1_FPOS$PES1 == "S 10") & (Beh1_FPOS$type == "S  9" | Beh1_FPOS$type == "S 10") & (Beh1_FPOS$prev == "S  2")), 0, NA))
Beh1_FPOS$PostError_con_EC <- ifelse(((Beh1_FPOS$PES1 == "S  9" | Beh1_FPOS$PES1 == "S 10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  1")), 1,
                                     ifelse(((Beh1_FPOS$PES1 == "S  9" | Beh1_FPOS$PES1 == "S 10") & (Beh1_FPOS$type == "S  9" | Beh1_FPOS$type == "S 10") & (Beh1_FPOS$prev == "S  1")), 0, NA))

#Pre Error accuracy
Beh1_FPOS$Pre_Error_inc_EC <- ifelse(((Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4") & (Beh1_FPOS$PESneg1 == "S  9" | Beh1_FPOS$PESneg1 == "S 10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  2")), 1,
                                     ifelse(((Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4") & (Beh1_FPOS$PESneg1 == "S  9" | Beh1_FPOS$PESneg1 == "S 10") & (Beh1_FPOS$type == "S  9" | Beh1_FPOS$type == "S 10") & (Beh1_FPOS$prev == "S  2")), 0, NA))
Beh1_FPOS$Pre_Error_con_EC <- ifelse(((Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4") & (Beh1_FPOS$PESneg1 == "S  9" | Beh1_FPOS$PESneg1 == "S 10") & (Beh1_FPOS$type == "S  3" | Beh1_FPOS$type == "S  4") & (Beh1_FPOS$prev == "S  1")), 1,
                                     ifelse(((Beh1_FPOS$PES1 == "S  3" | Beh1_FPOS$PES1 == "S  4") & (Beh1_FPOS$PESneg1 == "S  9" | Beh1_FPOS$PESneg1 == "S 10") & (Beh1_FPOS$type == "S  9" | Beh1_FPOS$type == "S 10") & (Beh1_FPOS$prev == "S  1")), 0, NA))

#Congruent number RTs
Beh1_FPOS$CongruentRT <- ifelse((Beh1_FPOS$prev == "S  1"), Beh1_FPOS$RT_trimmed, NA)

#Incongruent number RTs
Beh1_FPOS$IncongruentRT <- ifelse((Beh1_FPOS$prev == "S  2"), Beh1_FPOS$RT_trimmed, NA)

#Coding for Congruent vs. Incongruent and Correct vs. Error
Beh1_FPOS$Cong_Correct <- ifelse((Beh1_FPOS$prev == "S  1"), ifelse(Beh1_FPOS$EC == 1, 1, 0), NA)
Beh1_FPOS$Incong_Correct <- ifelse((Beh1_FPOS$prev == "S  2"), ifelse(Beh1_FPOS$EC == 1, 1, 0), NA)


#Calculate averages
BehAve_FPOS <- aggregate(Beh1_FPOS[, c("EC", "RT", "RT_trimmed", "CorrectRT", "ErrorRT", "CAfterCRT", "EAfterCRT", "CAfterERT", 
                                       "PostError_inc", "PostError_con", "Pre_Error_inc", "Pre_Error_con",
                                       "CongruentRT", "IncongruentRT", "Cong_Correct", "Incong_Correct")], list(Beh1_FPOS$ID), mean, na.rm=TRUE)

BehAve_FPOS$PES_Trad <- BehAve_FPOS$CAfterERT - BehAve_FPOS$CAfterCRT
BehAve_FPOS$PES_Robust <- (BehAve_FPOS$PostError_inc + BehAve_FPOS$PostError_con)/2 - (BehAve_FPOS$Pre_Error_inc + BehAve_FPOS$Pre_Error_con)/2


############################ 
############################ SUMMARY STATS & EXPORTING DATA
############################

print(describe(BehAve_FARR), digits = 4)
summary(BehAve_FARR)
hist(BehAve_FARR$RT)

print(describe(BehAve_FNEG), digits = 4)
summary(BehAve_FNEG)
hist(BehAve_FNEG$RT)

print(describe(BehAve_FPOS), digits = 4)
summary(BehAve_FPOS)
hist(BehAve_FPOS$RT)

## rename+ trim IDs to 4-digit # only so they match so can merge together & w/ EEG data later
BehAve_FARR$ID <- substr(BehAve_FARR$Group.1, 10, 13)
BehAve_FNEG$ID <- substr(BehAve_FNEG$Group.1, 10, 13)
BehAve_FPOS$ID <- substr(BehAve_FPOS$Group.1, 10, 13)

## export separately to then merge in with corresponding EEG data and select based on error_n 
write.csv(BehAve_FARR, file = paste("./Behav Extract Data/BehAve_FARR.csv"), na = "", row.names = FALSE)
write.csv(BehAve_FNEG, file = paste("./Behav Extract Data/BehAve_FNEG.csv"), na = "", row.names = FALSE)
write.csv(BehAve_FPOS, file = paste("./Behav Extract Data/BehAve_FPOS.csv"), na = "", row.names = FALSE)



