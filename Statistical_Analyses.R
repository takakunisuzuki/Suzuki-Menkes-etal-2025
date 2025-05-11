#####################################################################
# Response Monitoring Theta-Band Activities across Emotional Contexts in Schizophrenia- and Bipolar-Spectrum Disorders
# Suzuki, Menkes, et al.
#
# Written by Margo W. Menkes with help from Takakuni Suzuki
#####################################################################
rm(list=ls())

library(dplyr)
library(psych)
library(reshape2)
library(lme4)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(effectsize)

#################################################
############# DATA SETUP
###################################################
setwd("") # Set path

# Read in EEG data
FARR <- read.csv("./Output_Auto_Combined/EEG_FARR_Data_Resp_Combined_2025-05-11.csv", header=T, sep=",")
FNEG <- read.csv("./Output_Auto_Combined/EEG_FNEG_Data_Resp_Combined_2025-05-11.csv", header=T, sep=",")
FPOS <- read.csv("./Output_Auto_Combined/EEG_FPOS_Data_Resp_Combined_2025-05-11.csv", header=T, sep=",")

FARR_EBL <- read.csv("./Output_Auto_Combined/EEG_FARR_Data_Resp_EBL_2025-05-11.csv", header=T, sep=",")
FNEG_EBL <- read.csv("./Output_Auto_Combined/EEG_FNEG_Data_Resp_EBL_2025-05-11.csv", header=T, sep=",")
FPOS_EBL <- read.csv("./Output_Auto_Combined/EEG_FPOS_Data_Resp_EBL_2025-05-11.csv", header=T, sep=",")

# Read in behavioral data extracted/output from R
BehAve_FARR <- read.csv("./Behav Extract Data/BehAve_FARR.csv", header=T, sep=",")
BehAve_FNEG <- read.csv("./Behav Extract Data/BehAve_FNEG.csv", header=T, sep=",")
BehAve_FPOS <- read.csv("./Behav Extract Data/BehAve_FPOS.csv", header=T, sep=",")

# Read in participant characteristics and trim
subchar <- read.csv("PsyTheta_DATA_2025-03-27.csv", header=T, na.strings = c(999), sep=",")
subchar2 <- subchar[c("record_id", "demo_subjtype", "demo_age", "demo_sex", "demo_gender", "demo_race___1", "demo_race___2", "demo_race___3", "demo_race___4", "demo_race___5", "demo_race___6",
                      "demo_race___7", "demo_race___8", "demo_yrs_of_education", "demo_ethnicity", "demo_parent_education_yrs", "scid_primary_diagnosis", "scid_bpi_pf", "scid_psychosis_aao",
                      "hamd_totalscore", "hamd_totalscore_2nd",  "ymrs_totalscore", "ymrs_totalscore_2nd", "panss_positivescale_total", "panss_positivescale_total_2nd", 
                      "panss_negativescale_total", "panss_negativescale_total_2nd", "panss_generalpsychopathology_total", "panss_generalpsychopathology_total_2nd",
                      "hamd_repeat", "ymrs_repeat", "panss_repeat",
                      "demo_atyp_antipsy", "demo_typ_antipsy", "demo_anticonvuls", "demo_stim", "demo_bzd", "demo_ssri_snri", "demo_lith",
                      "farr_vas_p", "farr_vas_u", "farr_vas_e",
                      "fneg_vas_p", "fneg_vas_u", "fneg_vas_e", 
                      "fpos_vas_p", "fpos_vas_u", "fpos_vas_e", 
                      "scid_bpi_pf")]

# Merge EEG files
FARR <- merge(FARR, FARR_EBL[,c("ID", "ERN", "CRN", "ERN_non", "CRN_non")], by = c("ID"), suffixes = c("", ".ebl"), all = TRUE)
FNEG <- merge(FNEG, FNEG_EBL[,c("ID", "ERN", "CRN", "ERN_non", "CRN_non")], by = c("ID"), suffixes = c("", ".ebl"), all = TRUE)
FPOS <- merge(FPOS, FPOS_EBL[,c("ID", "ERN", "CRN", "ERN_non", "CRN_non")], by = c("ID"), suffixes = c("", ".ebl"), all = TRUE)

# Trim IDs to 4-digit # only so they match to merge w/ other data based on IDs
FARR$ID <- substr(FARR$ID, 10, 13)
FNEG$ID <- substr(FNEG$ID, 10, 13)
FPOS$ID <- substr(FPOS$ID, 10, 13)
subchar2$ID <- substr(subchar$record_id, 4, 7)

# Merge EEG and behavioral data
merged_FARR <- merge(FARR, BehAve_FARR, by = c("ID"), all = TRUE)
merged_FNEG <- merge(FNEG, BehAve_FNEG, by = c("ID"), all = TRUE)
merged_FPOS <- merge(FPOS, BehAve_FPOS, by = c("ID"), all = TRUE)

# Only keep those with at least 8 errors and 65% accuracy
FARR_wide <- subset(merged_FARR, Error_n >= 8 & Accuracy > .65)
FNEG_wide<- subset(merged_FNEG, Error_n >= 8 & Accuracy > .65)
FPOS_wide <- subset(merged_FPOS, Error_n >= 8 & Accuracy > .65)

# label variables within each data  by which task they are 
col_names_FARR <- names(FARR_wide)
exc_columns <- c("ID", "Group")
col_to_rename_FARR <- col_names_FARR[!(col_names_FARR %in% exc_columns)]
names(FARR_wide) <- ifelse(col_names_FARR %in% exc_columns, col_names_FARR, paste0("FARR_", col_names_FARR))

col_names_FNEG <- names(FNEG_wide)
col_to_rename_FNEG <- col_names_FNEG[!(col_names_FNEG %in% exc_columns)]
names(FNEG_wide) <- ifelse(col_names_FNEG %in% exc_columns, col_names_FNEG, paste0("FNEG_", col_names_FNEG))

col_names_FPOS <- names(FPOS_wide)
col_to_rename_FPOS <- col_names_FPOS[!(col_names_FPOS %in% exc_columns)]
names(FPOS_wide) <- ifelse(col_names_FPOS %in% exc_columns, col_names_FPOS, paste0("FPOS_", col_names_FPOS))

# merge EEG DFs
EEG_wide <- merge(FARR_wide, FNEG_wide, by = c("ID", "Group"), all = TRUE)
EEG_wide <- merge(EEG_wide, FPOS_wide, by = c("ID", "Group"), all = TRUE) 

# merge in key characteristics  
All_wide <- merge(EEG_wide, subchar2, by = c("ID")) 

#relabel group var also in larger df with full subj characteristics 
All_wide$Group[All_wide$Group==1] <- "HC"
All_wide$Group[All_wide$Group==2] <- "SZ"
All_wide$Group[All_wide$Group==3] <- "BD"
All_wide$Group <- factor(All_wide$Group, levels = c("HC", "BD", "SZ"))


#### Participant Characteristic Setup
### DEMO variable setup
# create labeled bio sex var
All_wide$sex_l <- NA
All_wide$sex_l[All_wide$demo_sex==0] <- "M"
All_wide$sex_l[All_wide$demo_sex==1] <- "F"

# create race var
All_wide$race_l <- NA
All_wide$race_l[All_wide$demo_race___1==1] <- "White"
All_wide$race_l[All_wide$demo_race___2==1] <- "B or AA"
All_wide$race_l[All_wide$demo_race___3==1] <- "Asian"
All_wide$race_l[All_wide$demo_race___4==1] <- "Native Am"
All_wide$race_l[All_wide$demo_race___5==1] <- "Native H or PI"
All_wide$race_l[All_wide$demo_race___6==1] <- "Multiple"
All_wide$race_l[All_wide$demo_race___7==1] <- "Unknown"
All_wide$race_l[All_wide$demo_race___8==1] <- "Other"

## collapse race into 2 grp (white vs. other/nonwhite race, bc counts are small in nonwhite groups) 
All_wide$Race_2grps[All_wide$demo_race___1==1] <- "White"
All_wide$Race_2grps[All_wide$demo_race___2==1|All_wide$demo_race___3==1|All_wide$demo_race___4==1|All_wide$demo_race___5==1|All_wide$demo_race___6==1|All_wide$demo_race___8==1] <- "Nonwhite"
All_wide$Race_2grps[All_wide$demo_race___7==1] <- "Unknown"

## create labeled ethnicity var
All_wide$ethnic_l <- NA
All_wide$ethnic_l[All_wide$demo_ethnicity==0] <- "Non H/L"
All_wide$ethnic_l[All_wide$demo_ethnicity==1] <- "H/L"

## create cleaned up gender identity var (this var was free text entry, inconsistencies in data entry)
All_wide$gender_cl <- NA
All_wide$gender_cl[All_wide$demo_gender=="female"] <- "F"
All_wide$gender_cl[All_wide$demo_gender=="male"|All_wide$demo_gender=="male (transition from F ->M)"|All_wide$demo_gender=="male (trans F ->M)"|All_wide$demo_gender=="trans male"] <- "M"
All_wide$gender_cl[All_wide$demo_gender=="nonbinary"|All_wide$demo_gender=="agender"|All_wide$demo_gender=="agender/nonbinary"|All_wide$demo_gender=="genderfluid, she/they"] <- "Other"

#### clinical variable setup

## created labeled diagnosis var
All_wide$Dx_l <- NA
All_wide$Dx_l[All_wide$scid_primary_diagnosis==1] <- "HC"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==2] <- "SZ"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==3] <- "SZFORM"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==4] <- "SZAFF"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==5] <- "DELUS"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==6] <- "BRIEFPSY"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==7] <- "PSYNOS"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==8] <- "BP1"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==9] <- "BP2"
All_wide$Dx_l[All_wide$scid_primary_diagnosis==10] <- "BPNOS"

## set up clin scales to use 2nd administration for any clin scales that were repeated clin scales (was repeated if 2wks+ passed since remote visit - want to use clin scales w/in 2 wks of EEG)
All_wide$hamd_use <- All_wide$hamd_totalscore
hamd_na_indices <- is.na(All_wide$hamd_totalscore_2nd) | is.na(All_wide$hamd_repeat) # Check for NAs in hamd_totalscore_2nd and hamd_repeat
All_wide$hamd_use[All_wide$hamd_repeat == 1 & !hamd_na_indices] <- All_wide$hamd_totalscore_2nd[!hamd_na_indices] # Update hamd_use only where hamd_repeat is 1 and hamd_totalscore_2nd is not NA

All_wide$ymrs_use <- All_wide$ymrs_totalscore
ymrs_na_indices <- is.na(All_wide$ymrs_totalscore_2nd) | is.na(All_wide$ymrs_repeat) # Check for NAs in ymrs_totalscore_2nd and ymrs_repeat
All_wide$ymrs_use[All_wide$ymrs_repeat == 1 & !ymrs_na_indices] <- All_wide$ymrs_totalscore_2nd[!ymrs_na_indices] # Update ymrs_use only where ymrs_repeat is 1 and ymrs_totalscore_2nd is not NA

All_wide$panss_pos_use <- All_wide$panss_positivescale_total
panss_pos_na_indices <- is.na(All_wide$panss_positivescale_total_2nd) | is.na(All_wide$panss_repeat) # Check for NAs in panss_positivescale_total_2nd and panss_repeat
All_wide$panss_pos_use[All_wide$panss_repeat == 1 & !panss_pos_na_indices] <- All_wide$panss_positivescale_total_2nd[!panss_pos_na_indices] # Update panss_pos_use only where panss_repeat is 1 and panss_positivescale_total_2nd is not NA

All_wide$panss_neg_use <- All_wide$panss_negativescale_total
panss_neg_na_indices <- is.na(All_wide$panss_negativescale_total_2nd) | is.na(All_wide$panss_repeat) # Check for NAs in panss_negativescale_total_2nd and panss_repeat
All_wide$panss_neg_use[All_wide$panss_repeat == 1 & !panss_neg_na_indices] <- All_wide$panss_negativescale_total_2nd[!panss_neg_na_indices] # Update panss_neg_use only where panss_repeat is 1 and panss_negativescale_total_2nd is not NA

All_wide$panss_gen_use <- All_wide$panss_generalpsychopathology_total
panss_gen_na_indices <- is.na(All_wide$panss_generalpsychopathology_total_2nd) | is.na(All_wide$panss_repeat) ## Check for NAs in panss_generalpsychopathology_total_2nd and panss_repeat
All_wide$panss_gen_use[All_wide$panss_repeat == 1 & !panss_gen_na_indices] <- All_wide$panss_generalpsychopathology_total_2nd[!panss_gen_na_indices] # Update panss_gen_use only where panss_repeat is 1 and panss_generalpsychopathology_total_2nd is not NA

## create var for mood state based on HAMD/YMRS
All_wide$mood_category <- NA
All_wide$mood_category[All_wide$hamd_use < 8 & All_wide$ymrs_use < 8] <- "euthym"
All_wide$mood_category[All_wide$hamd_use >= 8 & All_wide$ymrs_use < 8] <- "dep"
All_wide$mood_category[All_wide$hamd_use < 8 & All_wide$ymrs_use >=  8] <- "man"
All_wide$mood_category[All_wide$hamd_use >= 8 & All_wide$ymrs_use >=  8] <- "mixed"


# setup meds data so that 1 = yes (taking med), 2 = no (not taking) for all
# for each med type, 1 = yes (taking the med), 2 = no (not taking), blank if indicated 'no' to any psych meds in "demo_medication" var 
All_wide$any_antipsychot_use <- "2" # create variable combining typical & atypical antipsychotic to indicate if on any antipsychotic
All_wide$any_antipsychot_use [All_wide$demo_atyp_antipsy==1|All_wide$demo_typ_antipsy==1 ] <- "1"
All_wide$any_antipsychot_use [All_wide$demo_atyp_antipsy==2 & All_wide$demo_typ_antipsy==2 ] <- "2"

All_wide$demo_lith_use <- ifelse(is.na(All_wide$demo_lith), 2, All_wide$demo_lith) # create new lithium var that says they are not taking lithium if value is NA, otherwise keeps same as initial lithium var
All_wide$demo_anticonvuls_use <- ifelse(is.na(All_wide$demo_anticonvuls), 2, All_wide$demo_anticonvuls) # same for anticonvulsants 
All_wide$demo_bzd_use <- ifelse(is.na(All_wide$demo_bzd), 2, All_wide$demo_bzd) # same for benzos
All_wide$demo_stim_use <- ifelse(is.na(All_wide$demo_stim), 2, All_wide$demo_stim) # same for stimulants
All_wide$demo_ssri_snri_use <- ifelse(is.na(All_wide$demo_ssri_snri), 2, All_wide$demo_ssri_snri) # same for ssris/snris

# setup affect self-rating of each task
All_wide$farr_vas_p[All_wide$farr_vas_p==9999] <- NA
All_wide$farr_vas_u[All_wide$farr_vas_u==9999] <- NA
All_wide$farr_vas_e[All_wide$farr_vas_e==9999] <- NA
All_wide$fneg_vas_p[All_wide$fneg_vas_p==9999] <- NA
All_wide$fneg_vas_u[All_wide$fneg_vas_u==9999] <- NA
All_wide$fneg_vas_e[All_wide$fneg_vas_e==9999] <- NA
All_wide$fpos_vas_p[All_wide$fpos_vas_p==9999] <- NA
All_wide$fpos_vas_u[All_wide$fpos_vas_u==9999] <- NA
All_wide$fpos_vas_e[All_wide$fpos_vas_e==9999] <- NA

VAS <- reshape2::melt(All_wide[, c("ID", "Group", "farr_vas_p", "farr_vas_u", "farr_vas_e", "fneg_vas_p", "fneg_vas_u", "fneg_vas_e", "fpos_vas_p", "fpos_vas_u", "fpos_vas_e")], id.vars = c("ID", "Group"))
VAS<-VAS %>% dplyr::rename(rating = value)

VAS$task[VAS$variable=="farr_vas_p"|VAS$variable=="farr_vas_u"|VAS$variable=="farr_vas_e"] <- "FARR"
VAS$task[VAS$variable=="fneg_vas_p"|VAS$variable=="fneg_vas_u"|VAS$variable=="fneg_vas_e"] <- "FNEG"
VAS$task[VAS$variable=="fpos_vas_p"|VAS$variable=="fpos_vas_u"|VAS$variable=="fpos_vas_e"] <- "FPOS"

VAS$rating_type[VAS$variable=="farr_vas_p"|VAS$variable=="fneg_vas_p"|VAS$variable=="fpos_vas_p"] <- "pleas"
VAS$rating_type[VAS$variable=="farr_vas_u"|VAS$variable=="fneg_vas_u"|VAS$variable=="fpos_vas_u"] <- "unpleas"
VAS$rating_type[VAS$variable=="farr_vas_e"|VAS$variable=="fneg_vas_e"|VAS$variable=="fpos_vas_e"] <- "excit"

#### EEG data setup

# compute ERP/EEG diff scores within task to also show descriptives by group 
All_wide$FARR_delta_ERN <- (All_wide$FARR_ERN- All_wide$FARR_CRN)
All_wide$FARR_delta_ERN_non <- (All_wide$FARR_ERN_non - All_wide$FARR_CRN_non)
All_wide$FARR_delta_Theta_Power <- (All_wide$FARR_Power_Theta_E- All_wide$FARR_Power_Theta_C)
All_wide$FARR_delta_Theta_ITPC <- (All_wide$FARR_ITPC_Theta_E- All_wide$FARR_ITPC_Theta_C)

All_wide$FNEG_delta_ERN <- (All_wide$FNEG_ERN- All_wide$FNEG_CRN)
All_wide$FNEG_delta_ERN_non <- (All_wide$FNEG_ERN_non - All_wide$FNEG_CRN_non)
All_wide$FNEG_delta_Theta_Power <- (All_wide$FNEG_Power_Theta_E- All_wide$FNEG_Power_Theta_C)
All_wide$FNEG_delta_Theta_ITPC <- (All_wide$FNEG_ITPC_Theta_E- All_wide$FNEG_ITPC_Theta_C)

All_wide$FPOS_delta_ERN <- (All_wide$FPOS_ERN- All_wide$FPOS_CRN)
All_wide$FPOS_delta_ERN_non <- (All_wide$FPOS_ERN_non - All_wide$FPOS_CRN_non)
All_wide$FPOS_delta_Theta_Power <- (All_wide$FPOS_Power_Theta_E- All_wide$FPOS_Power_Theta_C)
All_wide$FPOS_delta_Theta_ITPC <- (All_wide$FPOS_ITPC_Theta_E- All_wide$FPOS_ITPC_Theta_C)


### Restructure into long format dataframes for each outcome variable 
# FARR
FARR_ERN <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ERN", "FARR_CRN", "FARR_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FARR_Accuracy"))
FARR_ERN_non <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ERN_non", "FARR_CRN_non", "FARR_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FARR_Accuracy"))
FARR_power <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FARR_Accuracy"))
FARR_itpc <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C", "FARR_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FARR_Accuracy"))
FARR_ERN_EBL <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ERN.ebl", "FARR_CRN.ebl", "FARR_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FARR_Accuracy"))

FARR_ERN<-FARR_ERN %>% dplyr::rename(accuracy = FARR_Accuracy, erp_amp = value, resp_type = variable)
FARR_ERN_non<-FARR_ERN_non %>% dplyr::rename(accuracy = FARR_Accuracy, erp_amp = value, resp_type = variable)
FARR_power<-FARR_power %>% dplyr::rename(accuracy = FARR_Accuracy, power = value, resp_type = variable)
FARR_itpc <-FARR_itpc %>% dplyr::rename(accuracy = FARR_Accuracy, itpc = value, resp_type = variable)
FARR_ERN_EBL<-FARR_ERN_EBL %>% dplyr::rename(accuracy = FARR_Accuracy, erp_amp = value, resp_type = variable)


# FNEG
FNEG_ERN<- reshape2::melt(All_wide[, c("ID", "Group", "FNEG_ERN", "FNEG_CRN", "FNEG_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FNEG_Accuracy"))
FNEG_ERN_non <- reshape2::melt(All_wide[, c("ID", "Group", "FNEG_ERN_non", "FNEG_CRN_non", "FNEG_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FNEG_Accuracy"))
FNEG_power <- reshape2::melt(All_wide[, c("ID", "Group", "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FNEG_Accuracy"))
FNEG_itpc <- reshape2::melt(All_wide[, c("ID", "Group", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C", "FNEG_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FNEG_Accuracy"))
FNEG_ERN_EBL<- reshape2::melt(All_wide[, c("ID", "Group", "FNEG_ERN.ebl", "FNEG_CRN.ebl", "FNEG_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FNEG_Accuracy"))

FNEG_ERN<-FNEG_ERN %>% dplyr::rename(accuracy = FNEG_Accuracy, erp_amp = value, resp_type = variable)
FNEG_ERN_non<-FNEG_ERN_non %>% dplyr::rename(accuracy = FNEG_Accuracy, erp_amp = value, resp_type = variable)
FNEG_power<-FNEG_power %>% dplyr::rename(accuracy = FNEG_Accuracy, power = value, resp_type = variable)
FNEG_itpc <-FNEG_itpc %>% dplyr::rename(accuracy = FNEG_Accuracy, itpc = value, resp_type = variable)

# FPOS
FPOS_ERN <- reshape2::melt(All_wide[, c("ID", "Group", "FPOS_ERN", "FPOS_CRN", "FPOS_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FPOS_Accuracy"))
FPOS_ERN_non <- reshape2::melt(All_wide[, c("ID", "Group", "FPOS_ERN_non", "FPOS_CRN_non", "FPOS_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FPOS_Accuracy"))
FPOS_power<- reshape2::melt(All_wide[, c("ID", "Group", "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FPOS_Accuracy"))
FPOS_itpc <- reshape2::melt(All_wide[, c("ID", "Group", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C", "FPOS_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FPOS_Accuracy"))
FPOS_ERN_EBL <- reshape2::melt(All_wide[, c("ID", "Group", "FPOS_ERN.ebl", "FPOS_CRN.ebl", "FPOS_Accuracy", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex", "FPOS_Accuracy"))

FPOS_ERN<-FPOS_ERN %>% dplyr::rename(accuracy = FPOS_Accuracy, erp_amp = value, resp_type = variable)
FPOS_ERN_non<-FPOS_ERN_non %>% dplyr::rename(accuracy = FPOS_Accuracy, erp_amp = value, resp_type = variable)
FPOS_power<-FPOS_power %>% dplyr::rename(accuracy = FPOS_Accuracy, power = value, resp_type = variable)
FPOS_itpc<-FPOS_itpc %>% dplyr::rename(accuracy = FPOS_Accuracy, itpc = value, resp_type = variable)
FPOS_ERN_EBL<-FPOS_ERN_EBL %>% dplyr::rename(accuracy = FPOS_Accuracy, erp_amp = value, resp_type = variable)


#### Restructure into long format dataframe with all data, including task valence as a factor
All_ERN <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ERN", "FARR_CRN", "FNEG_ERN", "FNEG_CRN", "FPOS_ERN", "FPOS_CRN","demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex"))
All_ERN_non <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ERN_non", "FARR_CRN_non", "FNEG_ERN_non", "FNEG_CRN_non", "FPOS_ERN_non", "FPOS_CRN_non", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex"))
All_power <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex"))
All_itpc <- reshape2::melt(All_wide[, c("ID", "Group",  "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex"))
All_ERN_EBL <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ERN.ebl", "FARR_CRN.ebl", "FNEG_ERN.ebl", "FNEG_CRN.ebl", "FPOS_ERN.ebl", "FPOS_CRN.ebl","demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex"))
All_ERN_non_EBL <- reshape2::melt(All_wide[, c("ID", "Group", "FARR_ERN_non.ebl", "FARR_CRN_non.ebl", "FNEG_ERN_non.ebl", "FNEG_CRN_non.ebl", "FPOS_ERN_non.ebl", "FPOS_CRN_non.ebl","demo_age", "demo_sex")], id.vars = c("ID", "Group", "demo_age", "demo_sex"))

All_ERN<-All_ERN %>% dplyr::rename(erp_amp = value)
All_ERN_non<-All_ERN_non %>% dplyr::rename(erp_amp = value)
All_power<-All_power %>% dplyr::rename(power = value)
All_itpc <-All_itpc %>% dplyr::rename(itpc = value)
All_ERN_EBL<-All_ERN_EBL %>% dplyr::rename(erp_amp = value)
All_ERN_non_EBL<-All_ERN_non_EBL %>% dplyr::rename(erp_amp = value)

### label within subj factors (response type + task valence) in long format data
# ERN matched
All_ERN$resp_type[All_ERN$variable=="FARR_ERN"|All_ERN$variable=="FNEG_ERN"|All_ERN$variable=="FPOS_ERN"] <- "ERN"
All_ERN$resp_type[All_ERN$variable=="FARR_CRN"|All_ERN$variable=="FNEG_CRN"|All_ERN$variable=="FPOS_CRN"] <- "CRN"
All_ERN$task_val[All_ERN$variable=="FARR_ERN"|All_ERN$variable=="FARR_CRN"] <- "FARR"
All_ERN$task_val[All_ERN$variable=="FNEG_ERN"|All_ERN$variable=="FNEG_CRN"] <- "FNEG"
All_ERN$task_val[All_ERN$variable=="FPOS_ERN"|All_ERN$variable=="FPOS_CRN"] <- "FPOS"

All_ERN_EBL$resp_type[All_ERN_EBL$variable=="FARR_ERN.ebl"|All_ERN_EBL$variable=="FNEG_ERN.ebl"|All_ERN_EBL$variable=="FPOS_ERN.ebl"] <- "ERN.ebl"
All_ERN_EBL$resp_type[All_ERN_EBL$variable=="FARR_CRN.ebl"|All_ERN_EBL$variable=="FNEG_CRN.ebl"|All_ERN_EBL$variable=="FPOS_CRN.ebl"] <- "CRN.ebl"
All_ERN_EBL$task_val[All_ERN_EBL$variable=="FARR_ERN.ebl"|All_ERN_EBL$variable=="FARR_CRN.ebl"] <- "FARR"
All_ERN_EBL$task_val[All_ERN_EBL$variable=="FNEG_ERN.ebl"|All_ERN_EBL$variable=="FNEG_CRN.ebl"] <- "FNEG"
All_ERN_EBL$task_val[All_ERN_EBL$variable=="FPOS_ERN.ebl"|All_ERN_EBL$variable=="FPOS_CRN.ebl"] <- "FPOS"

# ERN nonmatched
All_ERN_non$resp_type[All_ERN_non$variable=="FARR_ERN_non"|All_ERN_non$variable=="FNEG_ERN_non"|All_ERN_non$variable=="FPOS_ERN_non"] <- "ERN_non"
All_ERN_non$resp_type[All_ERN_non$variable=="FARR_CRN_non"|All_ERN_non$variable=="FNEG_CRN_non"|All_ERN_non$variable=="FPOS_CRN_non"] <- "CRN_non"
All_ERN_non$task_val[All_ERN_non$variable=="FARR_ERN_non"|All_ERN_non$variable=="FARR_CRN_non"] <- "FARR"
All_ERN_non$task_val[All_ERN_non$variable=="FNEG_ERN_non"|All_ERN_non$variable=="FNEG_CRN_non"] <- "FNEG"
All_ERN_non$task_val[All_ERN_non$variable=="FPOS_ERN_non"|All_ERN_non$variable=="FPOS_CRN_non"] <- "FPOS"

# ERN nonmatched using early baseline (EBL)
All_ERN_non_EBL$resp_type[All_ERN_non_EBL$variable=="FARR_ERN_non.ebl"|All_ERN_non_EBL$variable=="FNEG_ERN_non.ebl"|All_ERN_non_EBL$variable=="FPOS_ERN_non.ebl"] <- "ERN_non.ebl"
All_ERN_non_EBL$resp_type[All_ERN_non_EBL$variable=="FARR_CRN_non.ebl"|All_ERN_non_EBL$variable=="FNEG_CRN_non.ebl"|All_ERN_non_EBL$variable=="FPOS_CRN_non.ebl"] <- "CRN_non.ebl"
All_ERN_non_EBL$task_val[All_ERN_non_EBL$variable=="FARR_ERN_non.ebl"|All_ERN_non_EBL$variable=="FARR_CRN_non.ebl"] <- "FARR"
All_ERN_non_EBL$task_val[All_ERN_non_EBL$variable=="FNEG_ERN_non.ebl"|All_ERN_non_EBL$variable=="FNEG_CRN_non.ebl"] <- "FNEG"
All_ERN_non_EBL$task_val[All_ERN_non_EBL$variable=="FPOS_ERN_non.ebl"|All_ERN_non_EBL$variable=="FPOS_CRN_non.ebl"] <- "FPOS"

# theta power
All_power$resp_type[All_power$variable=="FARR_Power_Theta_E"|All_power$variable=="FNEG_Power_Theta_E"|All_power$variable=="FPOS_Power_Theta_E"] <- "Power_E"
All_power$resp_type[All_power$variable=="FARR_Power_Theta_C"|All_power$variable=="FNEG_Power_Theta_C"|All_power$variable=="FPOS_Power_Theta_C"] <- "Power_C"

All_power$task_val[All_power$variable=="FARR_Power_Theta_E"|All_power$variable=="FARR_Power_Theta_C"] <- "FARR"
All_power$task_val[All_power$variable=="FNEG_Power_Theta_E"|All_power$variable=="FNEG_Power_Theta_C"] <- "FNEG"
All_power$task_val[All_power$variable=="FPOS_Power_Theta_E"|All_power$variable=="FPOS_Power_Theta_C"] <- "FPOS"

# theta ITPC
All_itpc$resp_type[All_itpc$variable=="FARR_ITPC_Theta_E"|All_itpc$variable=="FNEG_ITPC_Theta_E"|All_itpc$variable=="FPOS_ITPC_Theta_E"] <- "ITPC_E"
All_itpc$resp_type[All_itpc$variable=="FARR_ITPC_Theta_C"|All_itpc$variable=="FNEG_ITPC_Theta_C"|All_itpc$variable=="FPOS_ITPC_Theta_C"] <- "ITPC_C"

All_itpc$task_val[All_itpc$variable=="FARR_ITPC_Theta_E"|All_itpc$variable=="FARR_ITPC_Theta_C"] <- "FARR"
All_itpc$task_val[All_itpc$variable=="FNEG_ITPC_Theta_E"|All_itpc$variable=="FNEG_ITPC_Theta_C"] <- "FNEG"
All_itpc$task_val[All_itpc$variable=="FPOS_ITPC_Theta_E"|All_itpc$variable=="FPOS_ITPC_Theta_C"] <- "FPOS"

### Restructure accuracy into long format + merge w/ EEG
All_acc <- reshape2::melt(All_wide[, c("ID", "FARR_Accuracy", "FNEG_Accuracy", "FPOS_Accuracy")], id.vars = c("ID"))
All_acc <-All_acc %>% dplyr::rename(accuracy = value)

All_acc$task_val[All_acc$variable=="FARR_Accuracy"] <- "FARR"
All_acc$task_val[All_acc$variable=="FNEG_Accuracy"] <- "FNEG"
All_acc$task_val[All_acc$variable=="FPOS_Accuracy"] <- "FPOS"

All_ERN <- merge(All_ERN, All_acc, by = c("ID", "task_val"), all = TRUE)
All_ERN_non <- merge(All_ERN_non, All_acc, by = c("ID", "task_val"), all = TRUE)
All_power <- merge(All_power, All_acc, by = c("ID", "task_val"), all = TRUE)
All_itpc <- merge(All_itpc, All_acc, by = c("ID", "task_val"), all = TRUE)
All_ERN_EBL <- merge(All_ERN_EBL, All_acc, by = c("ID", "task_val"), all = TRUE)
All_ERN_non_EBL <- merge(All_ERN_non_EBL, All_acc, by = c("ID", "task_val"), all = TRUE)

### Restructure into long format dataframe with all data, including task + antipsychotic med use as a factor
All_ERN_AntiPsychotic <- reshape2::melt(All_wide[, c("ID", "Group", "any_antipsychot_use", "FARR_ERN", "FARR_CRN", "FNEG_ERN", "FNEG_CRN", "FPOS_ERN", "FPOS_CRN","demo_age", "demo_sex")], id.vars = c("ID", "Group", "any_antipsychot_use", "demo_age", "demo_sex"))
All_power_AntiPsychotic <- reshape2::melt(All_wide[, c("ID", "Group", "any_antipsychot_use", "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "any_antipsychot_use", "demo_age", "demo_sex"))
All_itpc_AntiPsychotic <- reshape2::melt(All_wide[, c("ID", "Group", "any_antipsychot_use",  "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C", "demo_age", "demo_sex")], id.vars = c("ID", "Group", "any_antipsychot_use", "demo_age", "demo_sex"))

All_ERN_AntiPsychotic<-All_ERN_AntiPsychotic %>% dplyr::rename(erp_amp = value)
All_power_AntiPsychotic<-All_power_AntiPsychotic %>% dplyr::rename(power = value)
All_itpc_AntiPsychotic <-All_itpc_AntiPsychotic %>% dplyr::rename(itpc = value)

#### label within subj factors (response type + task valence) in long format data data w/ antipsychotic use var

# ERN 
All_ERN_AntiPsychotic$resp_type[All_ERN_AntiPsychotic$variable=="FARR_ERN"|All_ERN_AntiPsychotic$variable=="FNEG_ERN"|All_ERN_AntiPsychotic$variable=="FPOS_ERN"] <- "ERN"
All_ERN_AntiPsychotic$resp_type[All_ERN_AntiPsychotic$variable=="FARR_CRN"|All_ERN_AntiPsychotic$variable=="FNEG_CRN"|All_ERN_AntiPsychotic$variable=="FPOS_CRN"] <- "CRN"
All_ERN_AntiPsychotic$task_val[All_ERN_AntiPsychotic$variable=="FARR_ERN"|All_ERN_AntiPsychotic$variable=="FARR_CRN"] <- "FARR"
All_ERN_AntiPsychotic$task_val[All_ERN_AntiPsychotic$variable=="FNEG_ERN"|All_ERN_AntiPsychotic$variable=="FNEG_CRN"] <- "FNEG"
All_ERN_AntiPsychotic$task_val[All_ERN_AntiPsychotic$variable=="FPOS_ERN"|All_ERN_AntiPsychotic$variable=="FPOS_CRN"] <- "FPOS"

# theta power
All_power_AntiPsychotic$resp_type[All_power_AntiPsychotic$variable=="FARR_Power_Theta_E"|All_power_AntiPsychotic$variable=="FNEG_Power_Theta_E"|All_power_AntiPsychotic$variable=="FPOS_Power_Theta_E"] <- "Power_E"
All_power_AntiPsychotic$resp_type[All_power_AntiPsychotic$variable=="FARR_Power_Theta_C"|All_power_AntiPsychotic$variable=="FNEG_Power_Theta_C"|All_power_AntiPsychotic$variable=="FPOS_Power_Theta_C"] <- "Power_C"
All_power_AntiPsychotic$task_val[All_power_AntiPsychotic$variable=="FARR_Power_Theta_E"|All_power_AntiPsychotic$variable=="FARR_Power_Theta_C"] <- "FARR"
All_power_AntiPsychotic$task_val[All_power_AntiPsychotic$variable=="FNEG_Power_Theta_E"|All_power_AntiPsychotic$variable=="FNEG_Power_Theta_C"] <- "FNEG"
All_power_AntiPsychotic$task_val[All_power_AntiPsychotic$variable=="FPOS_Power_Theta_E"|All_power_AntiPsychotic$variable=="FPOS_Power_Theta_C"] <- "FPOS"

# theta ITPC
All_itpc_AntiPsychotic$resp_type[All_itpc_AntiPsychotic$variable=="FARR_ITPC_Theta_E"|All_itpc_AntiPsychotic$variable=="FNEG_ITPC_Theta_E"|All_itpc_AntiPsychotic$variable=="FPOS_ITPC_Theta_E"] <- "ITPC_E"
All_itpc_AntiPsychotic$resp_type[All_itpc_AntiPsychotic$variable=="FARR_ITPC_Theta_C"|All_itpc_AntiPsychotic$variable=="FNEG_ITPC_Theta_C"|All_itpc_AntiPsychotic$variable=="FPOS_ITPC_Theta_C"] <- "ITPC_C"
All_itpc_AntiPsychotic$task_val[All_itpc_AntiPsychotic$variable=="FARR_ITPC_Theta_E"|All_itpc_AntiPsychotic$variable=="FARR_ITPC_Theta_C"] <- "FARR"
All_itpc_AntiPsychotic$task_val[All_itpc_AntiPsychotic$variable=="FNEG_ITPC_Theta_E"|All_itpc_AntiPsychotic$variable=="FNEG_ITPC_Theta_C"] <- "FNEG"
All_itpc_AntiPsychotic$task_val[All_itpc_AntiPsychotic$variable=="FPOS_ITPC_Theta_E"|All_itpc_AntiPsychotic$variable=="FPOS_ITPC_Theta_C"] <- "FPOS"

# merge in accuracy data 
All_ERN_AntiPsychotic <- merge(All_ERN_AntiPsychotic, All_acc, by = c("ID", "task_val"), all = TRUE)
All_power_AntiPsychotic <- merge(All_power_AntiPsychotic, All_acc, by = c("ID", "task_val"), all = TRUE)
All_itpc_AntiPsychotic <- merge(All_itpc_AntiPsychotic, All_acc, by = c("ID", "task_val"), all = TRUE)

# subset into SZ group only for later analyses of antipsychotic effects

All_ERN_AntiPsychotic_SZ <- All_ERN_AntiPsychotic[All_ERN_AntiPsychotic$Group == "SZ",]
All_power_AntiPsychotic_SZ <- All_power_AntiPsychotic[All_power_AntiPsychotic$Group == "SZ",]
All_itpc_AntiPsychotic_SZ <- All_itpc_AntiPsychotic[All_itpc_AntiPsychotic$Group == "SZ",]

# subset into BD group only for later analyses of antipsychotic effects

All_ERN_AntiPsychotic_BD <- All_ERN_AntiPsychotic[All_ERN_AntiPsychotic$Group == "BD",]
All_power_AntiPsychotic_BD <- All_power_AntiPsychotic[All_power_AntiPsychotic$Group == "BD",]
All_itpc_AntiPsychotic_BD <- All_itpc_AntiPsychotic[All_itpc_AntiPsychotic$Group == "BD",]

#################################################
############# DESCRIPTIVE STATS
############# Tables S2, S3, S4, S11, S12 
###################################################

#### Table S2: Demo + clinical characteristics 

### Demographic

psych::describeBy(All_wide[,c("demo_age", "demo_yrs_of_education", "demo_parent_education_yrs")], All_wide$Group)
anova(aov(demo_age ~ Group, data=All_wide))
anova(aov(demo_yrs_of_education ~ Group, data=All_wide))
anova(aov(demo_parent_education_yrs ~ Group, data=All_wide))

# sex
table(All_wide$Group, All_wide$sex_l)
chisq.test(All_wide$Group, All_wide$sex_l)

# gender
table(All_wide$Group, All_wide$gender_cl)
chisq.test(All_wide$Group, All_wide$gender_cl)

# race
table(All_wide$Group, All_wide$race_l)
chisq.test(All_wide$Group, All_wide$race_l)

# chi sq test for race using 2-group race var due to small cell sizes in nonwhite grps 
table(All_wide$Group, All_wide$Race_2grps)
chisq.test(All_wide$Group, All_wide$Race_2grps)

# ethnicity 
table(All_wide$Group, All_wide$ethnic_l)
chisq.test(All_wide$Group, All_wide$ethnic_l)

### Clinical Characteristics

# diagnosis
table(All_wide$Group, All_wide$Dx_l)

# clin scales 
psych::describeBy(All_wide[, c("hamd_use", "ymrs_use", "panss_pos_use", "panss_neg_use", "panss_gen_use", "scid_psychosis_aao")], All_wide$Group)

anova(aov(hamd_use ~ Group, data=All_wide))
anova(aov(ymrs_use ~ Group, data=All_wide))
anova(aov(panss_pos_use ~ Group, data=All_wide))
anova(aov(panss_neg_use ~ Group, data=All_wide))
anova(aov(panss_gen_use ~ Group, data=All_wide))

## t.tests for post hoc contrasts on clin scales w/ no correction to get uncorrected p values 

# hamd
t.test(hamd_use ~ as.factor(Group), var.equal=TRUE,data = All_wide[All_wide$Group == "HC" | All_wide$Group == "BD", ])
t.test(hamd_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "SZ", ])
t.test(hamd_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ])

# ymrs
t.test(ymrs_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "BD", ])
t.test(ymrs_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "SZ", ])
t.test(ymrs_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ])

# panss-pos
t.test(panss_pos_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "BD", ])
t.test(panss_pos_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "SZ", ])
t.test(panss_pos_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ])

# panss-neg
t.test(panss_neg_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "BD", ])
t.test(panss_neg_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "SZ", ])
t.test(panss_neg_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ])

# panss-gen
t.test(panss_gen_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "BD", ])
t.test(panss_gen_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "HC" | All_wide$Group == "SZ", ])
t.test(panss_gen_use ~ as.factor(Group), var.equal=TRUE,data=All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ])


### medications - use subsetted data with BD/SZ grps only (no HCs)
# now 1 = yes (Taking the med), 2 = no (not taking med) 
table(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$any_antipsychot_use)
table(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_lith_use)
table(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_anticonvuls_use)
table(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_ssri_snri_use)
table(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_stim_use)
table(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_bzd_use)

## test for diffs betw 2 patient groups in med use 
chisq.test(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$any_antipsychot_use)
chisq.test(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_lith_use)
chisq.test(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_anticonvuls_use)
chisq.test(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_ssri_snri_use)
chisq.test(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_stim_use)
chisq.test(All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$Group, All_wide[All_wide$Group == "BD" | All_wide$Group == "SZ", ]$demo_bzd_use)

#### EEG descriptives by task & by group
psych::describeBy(All_wide[, c("FARR_ERN", "FARR_CRN", "FARR_CRN_non", "FARR_delta_ERN", "FARR_delta_ERN_non",
                                  "FNEG_ERN", "FNEG_CRN", "FNEG_CRN_non", "FNEG_delta_ERN", "FNEG_delta_ERN_non",
                                  "FPOS_ERN", "FPOS_CRN", "FPOS_CRN_non", "FPOS_delta_ERN", "FPOS_delta_ERN_non",
                                  
                                  "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_delta_Theta_Power", 
                                  "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_delta_Theta_Power", 
                                  "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_delta_Theta_Power",
                                  
                                  "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C", "FARR_delta_Theta_ITPC",
                                  "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C", "FNEG_delta_Theta_ITPC",
                                  "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C", "FPOS_delta_Theta_ITPC")], All_wide$Group)

#### Behavioral descriptives by task & group
psych::describeBy(All_wide[, c("FARR_EC", "FARR_Cong_Correct", "FARR_Incong_Correct", "FARR_RT_trimmed", "FARR_CongruentRT", "FARR_IncongruentRT", 
                               "FARR_CorrectRT", "FARR_ErrorRT", "FARR_CAfterCRT", "FARR_EAfterCRT", "FARR_CAfterERT", "FARR_PES_Trad", "FARR_PES_Robust", 
                               "FARR_Error_n", "FARR_Correct_n",
                               
                               "FNEG_EC", "FNEG_Cong_Correct", "FNEG_Incong_Correct", "FNEG_RT_trimmed", "FNEG_CongruentRT", "FNEG_IncongruentRT", 
                               "FNEG_CorrectRT", "FNEG_ErrorRT", "FNEG_CAfterCRT", "FNEG_EAfterCRT", "FNEG_CAfterERT", "FNEG_PES_Trad", "FNEG_PES_Robust", 
                               "FNEG_Error_n", "FNEG_Correct_n",
                               
                               
                               "FPOS_EC", "FPOS_Cong_Correct", "FPOS_Incong_Correct", "FPOS_RT_trimmed", "FPOS_CongruentRT", "FPOS_IncongruentRT", 
                               "FPOS_CorrectRT", "FPOS_ErrorRT", "FPOS_CAfterCRT", "FPOS_EAfterCRT", "FPOS_CAfterERT", "FPOS_PES_Trad", "FPOS_PES_Robust",
                               "FPOS_Error_n", "FPOS_Correct_n")], All_wide$Group)

### Mood state 
table(All_wide$Group, All_wide$mood_category)

All_wide %>% group_by(Group) %>% summarise(
  count_euthymic = sum(hamd_use < 8 & ymrs_use < 8, na.rm = TRUE),
  count_hamd_abv_8 = sum(hamd_use >= 8, na.rm = TRUE),
  count_ymrs_abv_8 = sum(ymrs_use >= 8, na.rm = TRUE),
  count_both_abv_8 = sum(hamd_use >= 8 & ymrs_use >= 8, na.rm = TRUE)
)

#### Affect-self rating descriptives by group
psych::describeBy(All_wide[,c("farr_vas_p", "farr_vas_u", "farr_vas_e",
                              "fneg_vas_p", "fneg_vas_u", "fneg_vas_e", 
                              "fpos_vas_p", "fpos_vas_u", "fpos_vas_e")], All_wide$Group)


#################################################
############# Inferential Statistics
###################################################

# Creating planned contrast list to test specific contrasts
contrast_list <- list(
  contrast_1 = c(2,-1,-1), # HC vs BD/SZ
  contrast_2 = c(0,1,-1) # BD vs SZ
)


## 2-way (Task x Response) ANCOVA within ea. group

# HC
HC_ERN_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_ERN[All_ERN$Group == "HC",])
anova(HC_ERN_lm)
emmeans(HC_ERN_lm, pairwise ~ task_val, adjust="none")
effectsize::eta_squared(HC_ERN_lm, partial = TRUE)

HC_power_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_power[All_power$Group == "HC",])
anova(HC_power_lm)
effectsize::eta_squared(HC_power_lm, partial = TRUE)

HC_itpc_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_itpc[All_itpc$Group == "HC",])
anova(HC_itpc_lm)
emmeans(HC_itpc_lm, pairwise ~ task_val | resp_type, adjust="none")
emmeans(HC_itpc_lm, pairwise ~ resp_type | task_val, adjust="none")
effectsize::eta_squared(HC_itpc_lm, partial = TRUE)

# BD
BD_ERN_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_ERN[All_ERN$Group == "BD",])
anova(BD_ERN_lm)
emmeans(BD_ERN_lm, pairwise ~ task_val, adjust="none")
effectsize::eta_squared(BD_ERN_lm, partial = TRUE)

BD_power_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_power[All_power$Group == "BD",])
anova(BD_power_lm)
effectsize::eta_squared(BD_power_lm, partial = TRUE)

BD_itpc_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_itpc[All_itpc$Group == "BD",])
anova(BD_itpc_lm)
emmeans(BD_itpc_lm, pairwise ~ task_val, adjust="none")
effectsize::eta_squared(BD_itpc_lm, partial = TRUE)

# SZ
SZ_ERN_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_ERN[All_ERN$Group == "SZ",])
anova(SZ_ERN_lm)
effectsize::eta_squared(SZ_ERN_lm, partial = TRUE)

SZ_power_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_power[All_power$Group == "SZ",])
anova(SZ_power_lm)
effectsize::eta_squared(SZ_power_lm, partial = TRUE)

SZ_itpc_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_itpc[All_itpc$Group == "SZ",])
anova(SZ_itpc_lm)
emmeans(SZ_itpc_lm, pairwise ~ task_val, adjust="none")
effectsize::eta_squared(SZ_itpc_lm, partial = TRUE)

## 2-way (Group x Response) ANCOVA within ea. task

# FARR
FARR_ERN_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy + resp_type*Group+(1|ID), REML = T, data =FARR_ERN)
anova(FARR_ERN_lm)
contrast(emmeans(FARR_ERN_lm, ~ Group, adjust="none"), contrast_list)
effectsize::eta_squared(FARR_ERN_lm, partial = TRUE)

FARR_power_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*Group+(1|ID), REML = T, data =FARR_power)
anova(FARR_power_lm)
contrast(emmeans(FARR_power_lm, ~ Group, adjust="none"), contrast_list)
emmeans(FARR_power_lm, pairwise ~ Group, adjust="none")
effectsize::eta_squared(FARR_power_lm, partial = TRUE)

FARR_itpc_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*Group+(1|ID), REML = T, data =FARR_itpc)
anova(FARR_itpc_lm)
contrast(emmeans(FARR_itpc_lm, ~ Group, adjust="none"), contrast_list)
effectsize::eta_squared(FARR_itpc_lm, partial = TRUE)

# FNEG
FNEG_ERN_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy + resp_type*Group+(1|ID), REML = T, data =FNEG_ERN)
anova(FNEG_ERN_lm)
contrast(emmeans(FNEG_ERN_lm, ~ Group, adjust="none"), contrast_list)
effectsize::eta_squared(FNEG_ERN_lm, partial = TRUE)

FNEG_power_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*Group+(1|ID), REML = T, data =FNEG_power)
anova(FNEG_power_lm)
contrast(emmeans(FNEG_power_lm, ~ Group, adjust="none"), contrast_list)
emmeans(FNEG_power_lm, pairwise ~ Group | resp_type, adjust="none")
emmeans(FNEG_power_lm, pairwise ~ resp_type | Group, adjust="none")
effectsize::eta_squared(FNEG_power_lm, partial = TRUE)

FNEG_itpc_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*Group+(1|ID), REML = T, data =FNEG_itpc)
anova(FNEG_itpc_lm)
contrast(emmeans(FNEG_itpc_lm, ~ Group, adjust="none"), contrast_list)
effectsize::eta_squared(FNEG_itpc_lm, partial = TRUE)

anova(lmerTest::lmer(itpc~demo_age+accuracy+resp_type+Group+(1|ID), REML = T, data =FNEG_itpc))
anova(lmerTest::lmer(itpc~demo_age+accuracy+Group+(1|ID), REML = T, data =FNEG_itpc))
anova(lmerTest::lmer(itpc~demo_age+accuracy+resp_type+(1|ID), REML = T, data =FNEG_itpc))
anova(lmerTest::lmer(itpc~demo_age+accuracy+(1|ID), REML = T, data =FNEG_itpc))
anova(lmerTest::lmer(itpc~demo_age+(1|ID), REML = T, data =FNEG_itpc))
anova(lmerTest::lmer(itpc~accuracy+(1|ID), REML = T, data =FNEG_itpc))
corr.test(FNEG_itpc$itpc, FNEG_itpc$accuracy)
anova(lmerTest::lmer(itpc~demo_age+resp_type*Group+(1|ID), REML = T, data =FNEG_itpc))
anova(lmerTest::lmer(itpc~demo_age+resp_type+Group+(1|ID), REML = T, data =FNEG_itpc))

# FPOS
FPOS_ERN_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy + resp_type*Group+(1|ID), REML = T, data =FPOS_ERN)
anova(FPOS_ERN_lm)
contrast(emmeans(FPOS_ERN_lm, ~ Group, adjust="none"), contrast_list)
effectsize::eta_squared(FPOS_ERN_lm, partial = TRUE)

FPOS_power_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*Group+(1|ID), REML = T, data =FPOS_power)
anova(FPOS_power_lm)
contrast(emmeans(FPOS_power_lm, ~ Group, adjust="none"), contrast_list)
emmeans(FPOS_power_lm, pairwise ~ Group | resp_type, adjust="none")
emmeans(FPOS_power_lm, pairwise ~ resp_type | Group, adjust="none")
effectsize::eta_squared(FPOS_power_lm, partial = TRUE)

FPOS_itpc_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*Group+(1|ID), REML = T, data =FPOS_itpc)
anova(FPOS_itpc_lm)
contrast(emmeans(FPOS_itpc_lm, ~ Group, adjust="none"), contrast_list)
effectsize::eta_squared(FPOS_itpc_lm, partial = TRUE)

anova(lmerTest::lmer(itpc~demo_age+accuracy+resp_type+Group+(1|ID), REML = T, data =FPOS_itpc))
anova(lmerTest::lmer(itpc~demo_age+accuracy+Group+(1|ID), REML = T, data =FPOS_itpc))
anova(lmerTest::lmer(itpc~demo_age+accuracy+resp_type+(1|ID), REML = T, data =FPOS_itpc))
anova(lmerTest::lmer(itpc~demo_age+accuracy+(1|ID), REML = T, data =FPOS_itpc))
anova(lmerTest::lmer(itpc~demo_age+(1|ID), REML = T, data =FPOS_itpc))
anova(lmerTest::lmer(itpc~accuracy+(1|ID), REML = T, data =FPOS_itpc))
corr.test(FPOS_itpc$itpc, FPOS_itpc$accuracy)
anova(lmerTest::lmer(itpc~demo_age+resp_type*Group+(1|ID), REML = T, data =FPOS_itpc))

#### 3-way 3X3X2 (Response × Group × Task) ANCOVA

## ERN 
All_ERN_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*Group*task_val+(1|ID), REML = T, data =All_ERN)
anova(All_ERN_lm)
contrast(emmeans(All_ERN_lm, ~ Group, adjust="none"), contrast_list)
emmeans(All_ERN_lm, pairwise ~ task_val, adjust="none")
effectsize::eta_squared(All_ERN_lm, partial = TRUE)

# store ERN results for later plotting of estimated marginal means
all_ERN_emm<-as.data.frame(emmeans(All_ERN_lm, ~ Group*resp_type*task_val, data=All_ERN))

## theta power
All_power_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*Group*task_val+(1|ID), REML = T, data =All_power)
anova(All_power_lm)
contrast(emmeans(All_power_lm, ~ Group, adjust="none"), contrast_list)
emmeans(All_power_lm, pairwise ~ Group | resp_type, adjust="none")
emmeans(All_power_lm, pairwise ~ resp_type | Group, adjust="none")
emmeans(All_power_lm, pairwise ~ task_val, adjust="none")
effectsize::eta_squared(All_power_lm, partial = TRUE)

# store theta power results for later plotting of estimated marginal means
all_power_emm<-as.data.frame(emmeans(All_power_lm, ~ resp_type*Group*task_val, data=All_power))

## theta ITPC
All_itpc_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*Group*task_val+(1|ID), REML = T, data =All_itpc)
anova(All_itpc_lm)
contrast(emmeans(All_itpc_lm, ~ Group, adjust="none"), contrast_list)
emmeans(All_itpc_lm, pairwise ~ task_val | resp_type, adjust="none")
emmeans(All_itpc_lm, pairwise ~ resp_type | task_val, adjust="none")
effectsize::eta_squared(All_itpc_lm, partial = TRUE)

# store theta ITPC results for later plotting of estimated marginal means
all_itpc_emm<-as.data.frame(emmeans(All_itpc_lm, ~ resp_type*Group*task_val, data=All_itpc))


#### Correlations

# within HC group correlations
corr_hc<-corr.test(All_wide[All_wide$Group == "HC", ][, c("FARR_ERN", "FARR_CRN", "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C", 
                             "FNEG_ERN", "FNEG_CRN", "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C",
                             "FPOS_ERN", "FPOS_CRN", "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C",                          
                             
                             "FARR_EC", "FARR_Cong_Correct", "FARR_Incong_Correct", "FARR_RT_trimmed", "FARR_CongruentRT", "FARR_IncongruentRT", 
                             "FARR_CorrectRT", "FARR_ErrorRT", "FARR_CAfterCRT", "FARR_EAfterCRT", "FARR_CAfterERT", "FARR_PES_Trad", "FARR_PES_Robust",

                             "FNEG_EC", "FNEG_Cong_Correct", "FNEG_Incong_Correct", "FNEG_RT_trimmed", "FNEG_CongruentRT", "FNEG_IncongruentRT", 
                             "FNEG_CorrectRT", "FNEG_ErrorRT", "FNEG_CAfterCRT", "FNEG_EAfterCRT", "FNEG_CAfterERT", "FNEG_PES_Trad", "FNEG_PES_Robust",
                             
                             "FPOS_EC", "FPOS_Cong_Correct", "FPOS_Incong_Correct", "FPOS_RT_trimmed", "FPOS_CongruentRT", "FPOS_IncongruentRT", 
                             "FPOS_CorrectRT", "FPOS_ErrorRT", "FPOS_CAfterCRT", "FPOS_EAfterCRT", "FPOS_CAfterERT", "FPOS_PES_Trad", "FPOS_PES_Robust")])  

write.csv(corr_hc$r, file = paste("./Correlation_1HC_checking.csv"), na = "", row.names = TRUE)

# within BD group correlations
corr_bd<-corr.test(All_wide[All_wide$Group == "BD", ][, c("FARR_ERN", "FARR_CRN", "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C", 
                             "FNEG_ERN", "FNEG_CRN", "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C",
                             "FPOS_ERN", "FPOS_CRN", "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C",                          
                             "FARR_EC", "FARR_Cong_Correct", "FARR_Incong_Correct", "FARR_RT_trimmed", "FARR_CongruentRT", "FARR_IncongruentRT", 
                             "FARR_CorrectRT", "FARR_ErrorRT", "FARR_CAfterCRT", "FARR_EAfterCRT", "FARR_CAfterERT", "FARR_PES_Trad", "FARR_PES_Robust",
                             "FNEG_EC", "FNEG_Cong_Correct", "FNEG_Incong_Correct", "FNEG_RT_trimmed", "FNEG_CongruentRT", "FNEG_IncongruentRT", 
                             "FNEG_CorrectRT", "FNEG_ErrorRT", "FNEG_CAfterCRT", "FNEG_EAfterCRT", "FNEG_CAfterERT", "FNEG_PES_Trad", "FNEG_PES_Robust",
                             "FPOS_EC", "FPOS_Cong_Correct", "FPOS_Incong_Correct", "FPOS_RT_trimmed", "FPOS_CongruentRT", "FPOS_IncongruentRT", 
                             "FPOS_CorrectRT", "FPOS_ErrorRT", "FPOS_CAfterCRT", "FPOS_EAfterCRT", "FPOS_CAfterERT", "FPOS_PES_Trad", "FPOS_PES_Robust",
                             "hamd_use", "ymrs_use", "panss_pos_use", "panss_neg_use", "panss_gen_use")])  

write.csv(corr_bd$r, file = paste("./Correlation_1BD_checking.csv"), na = "", row.names = TRUE)

# within SZ group correlations
corr_sz<-corr.test(All_wide[All_wide$Group == "SZ", ][, c("FARR_ERN", "FARR_CRN", "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C", 
                             "FNEG_ERN", "FNEG_CRN", "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C",
                             "FPOS_ERN", "FPOS_CRN", "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C",                          
                             
                             "FARR_EC", "FARR_Cong_Correct", "FARR_Incong_Correct", "FARR_RT_trimmed", "FARR_CongruentRT", "FARR_IncongruentRT", 
                             "FARR_CorrectRT", "FARR_ErrorRT", "FARR_CAfterCRT", "FARR_EAfterCRT", "FARR_CAfterERT", "FARR_PES_Trad", "FARR_PES_Robust",
                             "FNEG_EC", "FNEG_Cong_Correct", "FNEG_Incong_Correct", "FNEG_RT_trimmed", "FNEG_CongruentRT", "FNEG_IncongruentRT", 
                             "FNEG_CorrectRT", "FNEG_ErrorRT", "FNEG_CAfterCRT", "FNEG_EAfterCRT", "FNEG_CAfterERT", "FNEG_PES_Trad", "FNEG_PES_Robust",
                             "FPOS_EC", "FPOS_Cong_Correct", "FPOS_Incong_Correct", "FPOS_RT_trimmed", "FPOS_CongruentRT", "FPOS_IncongruentRT", 
                             "FPOS_CorrectRT", "FPOS_ErrorRT", "FPOS_CAfterCRT", "FPOS_EAfterCRT", "FPOS_CAfterERT", "FPOS_PES_Trad", "FPOS_PES_Robust",
                             
                             "hamd_use", "ymrs_use", "panss_pos_use", "panss_neg_use", "panss_gen_use")])  

write.csv(corr_sz$r, file = paste("./Correlation_1SZ_checking.csv"), na = "", row.names = TRUE)


#### antipsychotic effects on EEG

### antipsychotic effects on EEG within SZ
# ERN 
All_ERN_AntiPsychotic_SZ_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*any_antipsychot_use*task_val+(1|ID), REML = T, data =All_ERN_AntiPsychotic_SZ)
anova(All_ERN_AntiPsychotic_SZ_lm)

# theta power
All_power_AntiPsychotic_SZ_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*any_antipsychot_use*task_val+(1|ID), REML = T, data =All_power_AntiPsychotic_SZ)
anova(All_power_AntiPsychotic_SZ_lm)

# theta itpc
All_itpc_AntiPsychotic_SZ_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*any_antipsychot_use*task_val+(1|ID), REML = T, data =All_itpc_AntiPsychotic_SZ)
anova(All_itpc_AntiPsychotic_SZ_lm)

### antipsychotic effects on EEG within BD
# ERN  
All_ERN_AntiPsychotic_BD_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*any_antipsychot_use*task_val+(1|ID), REML = T, data =All_ERN_AntiPsychotic_BD)
anova(All_ERN_AntiPsychotic_BD_lm)

### theta power
All_power_AntiPsychotic_BD_lm<-lmerTest::lmer(power~demo_age+accuracy+resp_type*any_antipsychot_use*task_val+(1|ID), REML = T, data =All_power_AntiPsychotic_BD)
anova(All_power_AntiPsychotic_BD_lm)

## theta itpc
All_itpc_AntiPsychotic_BD_lm<-lmerTest::lmer(itpc~demo_age+accuracy+resp_type*any_antipsychot_use*task_val+(1|ID), REML = T, data =All_itpc_AntiPsychotic_BD)
anova(All_itpc_AntiPsychotic_BD_lm)


#### 1-way ANCOVA comparing affect self-rating for each task within each group

# within HC
VAS_HC_unpleas_lm<-aov(rating ~ task, data=VAS[VAS$Group == "HC" & VAS$rating_type == "unpleas", ])
anova(VAS_HC_unpleas_lm)
emmeans(VAS_HC_unpleas_lm, pairwise ~ task, adjust = "none")

VAS_HC_pleas_lm<-aov(rating ~ task, data=VAS[VAS$Group == "HC" & VAS$rating_type == "pleas", ])
anova(VAS_HC_pleas_lm)
emmeans(VAS_HC_pleas_lm, pairwise ~ task, adjust = "none")

VAS_HC_excit_lm<-aov(rating ~ task, data=VAS[VAS$Group == "HC" & VAS$rating_type == "excit", ])
anova(VAS_HC_excit_lm)

# within BD
VAS_BD_unpleas_lm<-aov(rating ~ task, data=VAS[VAS$Group == "BD" & VAS$rating_type == "unpleas", ])
anova(VAS_BD_unpleas_lm)
emmeans(VAS_BD_unpleas_lm, pairwise ~ task, adjust = "none")

VAS_BD_pleas_lm<-aov(rating ~ task, data=VAS[VAS$Group == "BD" & VAS$rating_type == "pleas", ])
anova(VAS_BD_pleas_lm)
emmeans(VAS_BD_pleas_lm, pairwise ~ task, adjust = "none")

VAS_BD_excit_lm<-aov(rating ~ task, data=VAS[VAS$Group == "BD" & VAS$rating_type == "excit", ])
anova(VAS_BD_excit_lm)

# within SZ
VAS_SZ_unpleas_lm<-aov(rating ~ task, data=VAS[VAS$Group == "SZ" & VAS$rating_type == "unpleas", ])
anova(VAS_SZ_unpleas_lm)
emmeans(VAS_SZ_unpleas_lm, pairwise ~ task, adjust = "none")

VAS_SZ_pleas_lm<-aov(rating ~ task, data=VAS[VAS$Group == "SZ" & VAS$rating_type == "pleas", ])
anova(VAS_SZ_pleas_lm)
emmeans(VAS_SZ_pleas_lm, pairwise ~ task, adjust = "none")

VAS_BD_excit_lm<-aov(rating ~ task, data=VAS[VAS$Group == "SZ" & VAS$rating_type == "excit", ])
anova(VAS_BD_excit_lm)

#### ANCOVAs of early-baseline ERPs 

## 3-way (Response × Group × Task) ANCOVA for early-baseline ERPs
All_ERN_EBL_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*Group*task_val+(1|ID), REML = T, data =All_ERN_non_EBL)
anova(All_ERN_EBL_non_lm)

## 2-way (Response X Task) ANCOVA for early-baseline ERPs within HC group.
HC_ERN_EBL_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_ERN_non_EBL[All_ERN_non_EBL$Group == "HC",])
anova(HC_ERN_EBL_non_lm)

#### 2-way ANCOVAs of ERPs incl. all trials

# Response X Group ANCOVA (within each task) for ERP incl. all trials

FARR_ERN_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy + resp_type*Group+(1|ID), REML = T, data =FARR_ERN_non)
anova(FARR_ERN_non_lm)
contrast(emmeans(FARR_ERN_non_lm, ~ Group, adjust="none"), contrast_list)

FNEG_ERN_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy + resp_type*Group+(1|ID), REML = T, data =FNEG_ERN_non)
anova(FNEG_ERN_non_lm)
contrast(emmeans(FNEG_ERN_non_lm, ~ Group, adjust="none"), contrast_list)

FPOS_ERN_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy + resp_type*Group+(1|ID), REML = T, data =FPOS_ERN_non)
anova(FPOS_ERN_non_lm)
contrast(emmeans(FPOS_ERN_non_lm, ~ Group, adjust="none"), contrast_list)

# Response X Task ANCOVAs (within each group) for ERP incl. all trials

HC_ERN_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_ERN_non[All_ERN_non$Group == "HC",])
anova(HC_ERN_non_lm)
emmeans(HC_ERN_non_lm, pairwise ~ task_val, adjust="none")

BD_ERN_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_ERN_non[All_ERN_non$Group == "BD",])
anova(BD_ERN_non_lm)
emmeans(BD_ERN_non_lm, pairwise ~ task_val, adjust="none")

SZ_ERN_non_lm<-lmerTest::lmer(erp_amp~demo_age+accuracy+resp_type*task_val+(1|ID), REML = T, data = All_ERN_non[All_ERN_non$Group == "SZ",])
anova(SZ_ERN_non_lm)


#################################################
############# Plotting
###################################################

#### rename for axes labels 

all_ERN_emm$Response <- recode_factor(all_ERN_emm$resp_type, ERN = "Error", CRN = "Correct")
all_ERN_emm$Task <- recode_factor(all_ERN_emm$task_val, FARR = "Arrow", FNEG = "Unpleasant", FPOS = "Pleasant")

all_power_emm$Response <- recode_factor(all_power_emm$resp_type, Power_E = "Error", Power_C = "Correct")
all_power_emm$Task <- recode_factor(all_power_emm$task_val, FARR = "Arrow", FNEG = "Unpleasant", FPOS = "Pleasant")

all_itpc_emm$Response <- recode_factor(all_itpc_emm$resp_type, ITPC_E = "Error", ITPC_C = "Correct")
all_itpc_emm$Task <- recode_factor(all_itpc_emm$task_val, FARR = "Arrow", FNEG = "Unpleasant", FPOS = "Pleasant")

#### Line plot of estimated marginal means from 3-way ANCOVA results

ERN_emms_lineplot <- ggplot(all_ERN_emm, aes(x = Task, y = emmean, color = Response, group = Response, linetype = Response)) +
  geom_line(linewidth = 1)  +  # lines connecting pts
  geom_point(position = position_dodge(width = 0.1)) +  # pts at each value
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(width = 0.1)) +  # error bars
  scale_color_manual(values = c("#999999", "#0066FF")) +
  labs(x = "", y = "ERP amplitude", fill = "Response", title = "ERP") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, margin = margin(t = 10)), ## adding space at bottom of plot (above x-axis text) to add in asterisks in ppt 
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))) +
  facet_wrap(~ Group)

power_emms_lineplot <- ggplot(all_power_emm, aes(x = Task, y = emmean, color = Response, group = Response, linetype = Response)) +
  geom_line(linewidth = 1) +  # lines connecting pts
  geom_point(position = position_dodge(width = 0.1)) +  # pts at each value
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(width = 0.1)) +  # error bars
  scale_color_manual(values = c("#999999", "#0066FF")) +
  labs(x = "", y = "Power (dB)", fill = "Response", title = "Power") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, margin = margin(t = 10)), ## adding space at bottom of plot (above x-axis text) to add in asterisks in ppt 
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))) +
  facet_wrap(~ Group)

itpc_emms_lineplot <- ggplot(all_itpc_emm, aes(x = Task, y = emmean, color = Response, group = Response, linetype = Response)) +
  geom_line(linewidth = 1) +  # lines connecting pts
  geom_point(position = position_dodge(width = 0.1)) +  # pts at each value
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(width = 0.1)) +  # error bars
  scale_color_manual(values = c("#999999", "#0066FF")) +
  labs(x = "", y = "ITPC", fill = "Response", title = "ITPC") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, margin = margin(t = 10)), ## adding space at bottom of plot (above x-axis text) to add in asterisks in ppt 
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        strip.text = element_text(size = 10, face = "bold", margin = margin(b = 10))) +
  facet_wrap(~ Group)

ggarrange(ERN_emms_lineplot, power_emms_lineplot, itpc_emms_lineplot,
          ncol = 1, nrow = 3, common.legend=TRUE, legend = "top")

#### Corr plots 

## ERN
Scatt_ERN_FARR_FNEG<-ggplot2::ggplot(All_wide, aes(x = FARR_ERN, y = FNEG_ERN, color = Group, shape = Group, linetype = Group)) +
  geom_point(size = 1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Arrow ERN", y = "Unpleasant ERN") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

Scatt_ERN_FARR_FPOS<-ggplot2::ggplot(All_wide, aes(x = FARR_ERN, y = FPOS_ERN, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Arrow ERN", y = "Pleasant ERN") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

Scatt_ERN_FNEG_FPOS<-ggplot2::ggplot(All_wide, aes(x = FNEG_ERN, y = FPOS_ERN, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Unpleasant ERN", y = "Pleasant ERN") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

## Error theta power
Scatt_theta_E_FARR_FNEG<-ggplot2::ggplot(All_wide, aes(x = FARR_Power_Theta_E, y = FNEG_Power_Theta_E, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Arrow Error Theta Power", y = "Unpleasant Error Theta Power") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

Scatt_theta_E_FARR_FPOS<-ggplot2::ggplot(All_wide, aes(x = FARR_Power_Theta_E, y = FPOS_Power_Theta_E, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +  # Add reg line
  labs(x = "Arrow Error Theta Power", y = "Pleasant Error Theta Power") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

Scatt_theta_E_FNEG_FPOS<-ggplot2::ggplot(All_wide, aes(x = FNEG_Power_Theta_E, y = FPOS_Power_Theta_E, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Unpleasant Error Theta Power", y = "Pleasant Error Theta Power") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

## Error Theta ITPC 
Scatt_ITPC_E_FARR_FNEG<-ggplot2::ggplot(All_wide, aes(x = FARR_ITPC_Theta_E, y = FNEG_ITPC_Theta_E, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Arrow Error Theta ITPC", y = "Unpleasant Error Theta ITPC") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

Scatt_ITPC_E_FARR_FPOS<-ggplot2::ggplot(All_wide, aes(x = FARR_ITPC_Theta_E, y = FPOS_ITPC_Theta_E, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Arrow Error Theta ITPC", y = "Pleasant Error Theta ITPC") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

Scatt_ITPC_E_FNEG_FPOS<-ggplot2::ggplot(All_wide, aes(x = FNEG_ITPC_Theta_E, y = FPOS_ITPC_Theta_E, color = Group, shape = Group, linetype = Group)) +
  geom_point(size=1.5) +  # Scatterplot
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # Add reg line
  labs(x = "Unpleasant Error Theta ITPC", y = "Pleasant Error Theta ITPC") + # Axis labels
  theme_minimal() + theme( axis.text = element_text(size = 10),  
                           axis.title = element_text(size = 10))+
  scale_color_manual(values = c("HC" = "#999999", "BD" = "#3366CC", "SZ" = "#486856"))

# create matrix of corr plots 
ggarrange(Scatt_ERN_FARR_FNEG, Scatt_ERN_FARR_FPOS, Scatt_ERN_FNEG_FPOS,
          Scatt_theta_E_FARR_FNEG, Scatt_theta_E_FARR_FPOS, Scatt_theta_E_FNEG_FPOS,
          Scatt_ITPC_E_FARR_FNEG, Scatt_ITPC_E_FARR_FPOS, Scatt_ITPC_E_FNEG_FPOS,
          ncol = 3, nrow = 3, common.legend=TRUE, legend = "top")


########
#Check <- merge(FARR_wide[,c("ID", "Accuracy")], FNEG_wide[,c("ID", "Accuracy")], by = "ID",suffixes = c(".arrow",".negative"), all = TRUE)
#Check <- merge(Check, FPOS_wide[,c("ID", "Accuracy")], by = "ID",suffixes = c("",".positive"), all = TRUE)
#write.csv(Check, file = "./check.csv", na = "")


####  SME 
FARR_SME <- read.csv("./Output_Auto_Combined/SME_FARR_2025-05-11.csv", header=T, sep=",")
FNEG_SME <- read.csv("./Output_Auto_Combined/SME_FNEG_2025-05-11.csv", header=T, sep=",")
FPOS_SME <- read.csv("./Output_Auto_Combined/SME_FPOS_2025-05-11.csv", header=T, sep=",")

FARR_SME <- subset(FARR_SME, Error_n >= 8)
FNEG_SME<- subset(FNEG_SME, Error_n >= 8)
FPOS_SME <- subset(FPOS_SME, Error_n >= 8)

FARR_SME$ID <- substr(FARR_SME$ID, 10, 13)
FNEG_SME$ID <- substr(FNEG_SME$ID, 10, 13)
FPOS_SME$ID <- substr(FPOS_SME$ID, 10, 13)

FARR_SNR <- merge(FARR_wide, FARR_SME, by = "ID")
FNEG_SNR <- merge(FNEG_wide, FNEG_SME, by = "ID")
FPOS_SNR <- merge(FPOS_wide, FPOS_SME, by = "ID")


# Calculating RMS
length(FARR_SNR$SME_theta_power_error)
sqrt(sum(FARR_SNR$SME_theta_power_error * FARR_SNR$SME_theta_power_error)/length(FARR_SNR$SME_theta_power_error))
sqrt(sum(FARR_SNR$SME_theta_power_correct * FARR_SNR$SME_theta_power_correct)/length(FARR_SNR$SME_theta_power_correct))
sqrt(sum(FARR_SNR$SME_ERN * FARR_SNR$SME_ERN)/length(FARR_SNR$SME_ERN))
sqrt(sum(FARR_SNR$SME_CRN * FARR_SNR$SME_CRN)/length(FARR_SNR$SME_CRN))
sqrt(sum(FARR_SNR$SME_CRN_non * FARR_SNR$SME_CRN_non)/length(FARR_SNR$SME_CRN_non))

length(FNEG_SNR$SME_theta_power_error)
sqrt(sum(FNEG_SNR$SME_theta_power_error * FNEG_SNR$SME_theta_power_error)/length(FNEG_SNR$SME_theta_power_error))
sqrt(sum(FNEG_SNR$SME_theta_power_correct * FNEG_SNR$SME_theta_power_correct)/length(FNEG_SNR$SME_theta_power_correct))
sqrt(sum(FNEG_SNR$SME_ERN * FNEG_SNR$SME_ERN)/length(FNEG_SNR$SME_ERN))
sqrt(sum(FNEG_SNR$SME_CRN * FNEG_SNR$SME_CRN)/length(FNEG_SNR$SME_CRN))
sqrt(sum(FNEG_SNR$SME_CRN_non * FNEG_SNR$SME_CRN_non)/length(FNEG_SNR$SME_CRN_non))

length(FPOS_SNR$SME_theta_power_error)
sqrt(sum(FPOS_SNR$SME_theta_power_error * FPOS_SNR$SME_theta_power_error)/length(FPOS_SNR$SME_theta_power_error))
sqrt(sum(FPOS_SNR$SME_theta_power_correct * FPOS_SNR$SME_theta_power_correct)/length(FPOS_SNR$SME_theta_power_correct))
sqrt(sum(FPOS_SNR$SME_ERN * FPOS_SNR$SME_ERN)/length(FPOS_SNR$SME_ERN))
sqrt(sum(FPOS_SNR$SME_CRN * FPOS_SNR$SME_CRN)/length(FPOS_SNR$SME_CRN))
sqrt(sum(FPOS_SNR$SME_CRN_non * FPOS_SNR$SME_CRN_non)/length(FPOS_SNR$SME_CRN_non))


# Calculate SNR
FARR_SNR$FARR_Power_Theta_E_SNR <- FARR_SNR$FARR_Power_Theta_E / FARR_SNR$SME_theta_power_error
FARR_SNR$FARR_ERN_SNR           <- FARR_SNR$FARR_ERN           / FARR_SNR$SME_ERN
FARR_SNR$FARR_Power_Theta_C_SNR <- FARR_SNR$FARR_Power_Theta_C / FARR_SNR$SME_theta_power_correct
FARR_SNR$FARR_CRN_SNR           <- FARR_SNR$FARR_CRN           / FARR_SNR$SME_CRN
FARR_SNR$FARR_CRN_non_SNR       <- FARR_SNR$FARR_CRN_non       / FARR_SNR$SME_CRN_non

FNEG_SNR$FNEG_Power_Theta_E_SNR <- FNEG_SNR$FNEG_Power_Theta_E / FNEG_SNR$SME_theta_power_error
FNEG_SNR$FNEG_ERN_SNR           <- FNEG_SNR$FNEG_ERN           / FNEG_SNR$SME_ERN
FNEG_SNR$FNEG_Power_Theta_C_SNR <- FNEG_SNR$FNEG_Power_Theta_C / FNEG_SNR$SME_theta_power_correct
FNEG_SNR$FNEG_CRN_SNR           <- FNEG_SNR$FNEG_CRN           / FNEG_SNR$SME_CRN
FNEG_SNR$FNEG_CRN_non_SNR       <- FNEG_SNR$FNEG_CRN_non       / FNEG_SNR$SME_CRN_non

FPOS_SNR$FPOS_Power_Theta_E_SNR <- FPOS_SNR$FPOS_Power_Theta_E / FPOS_SNR$SME_theta_power_error
FPOS_SNR$FPOS_ERN_SNR           <- FPOS_SNR$FPOS_ERN           / FPOS_SNR$SME_ERN
FPOS_SNR$FPOS_Power_Theta_C_SNR <- FPOS_SNR$FPOS_Power_Theta_C / FPOS_SNR$SME_theta_power_correct
FPOS_SNR$FPOS_CRN_SNR           <- FPOS_SNR$FPOS_CRN           / FPOS_SNR$SME_CRN
FPOS_SNR$FPOS_CRN_non_SNR       <- FPOS_SNR$FPOS_CRN_non       / FPOS_SNR$SME_CRN_non

SNR_all <- merge(FARR_SNR[,c("ID", "FARR_Power_Theta_E_SNR", "FARR_Power_Theta_C_SNR", "FARR_ERN_SNR", "FARR_CRN_SNR", "FARR_CRN_non_SNR")], 
                 FNEG_SNR[,c("ID", "FNEG_Power_Theta_E_SNR", "FNEG_Power_Theta_C_SNR", "FNEG_ERN_SNR", "FNEG_CRN_SNR", "FNEG_CRN_non_SNR")], 
                 by = "ID", all = TRUE)
SNR_all <- merge(SNR_all, FPOS_SNR[,c("ID", "FPOS_Power_Theta_E_SNR", "FPOS_Power_Theta_C_SNR", "FPOS_ERN_SNR", "FPOS_CRN_SNR", "FPOS_CRN_non_SNR")], 
                 by = "ID", all = TRUE)

describe(SNR_all[,-1])


### Cluster Analysis
All_wide_cluster <- All_wide[,c("ID", "Group", "scid_bpi_pf",
                                "FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C",
                                "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C",
                                "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C")]
All_wide_cluster[,c("FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C",
                    "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C",
                    "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C")] <- 
  scale(All_wide_cluster[,c("FARR_Power_Theta_E", "FARR_Power_Theta_C", "FARR_ITPC_Theta_E", "FARR_ITPC_Theta_C",
  "FNEG_Power_Theta_E", "FNEG_Power_Theta_C", "FNEG_ITPC_Theta_E", "FNEG_ITPC_Theta_C",
  "FPOS_Power_Theta_E", "FPOS_Power_Theta_C", "FPOS_ITPC_Theta_E", "FPOS_ITPC_Theta_C")])


All_wide_cluster$scid_bpi_pf <- ifelse(All_wide_cluster$Group == "SZ", 2,
                                       ifelse(All_wide_cluster$Group == "HC", 1, All_wide_cluster$scid_bpi_pf))

All_wide_cluster$scid_bpi_pf <- All_wide_cluster$scid_bpi_pf-1

# Listwise deletion
All_wide_cluster <- na.omit(All_wide_cluster)

Cluster_pc <- prcomp(All_wide_cluster[,-c(1:3)], scale = TRUE)
summary(Cluster_pc)
Cluster_pc$rotation[,1:3]

parallel = fa.parallel(All_wide_cluster[,-c(1:3)], fm = 'ml', fa = 'pc', nfactors = 3, n.iter = 100)

set.seed(1234)
n <- nrow(All_wide_cluster[,-c(1:3)])
Cluster_group <- matrix(NA, nrow = n, ncol = 3)
for(i in 2:3) {
  Cluster_group[,i] <- (kmeans(All_wide_cluster[,-c(1:3)], centers = i)$cluster)
}

Cluster_group[,1] <- All_wide_cluster$ID
Cluster_group <- merge(Cluster_group, All_wide_cluster[,c("ID", "Group", "scid_bpi_pf")], by.x = "V1", by.y = "ID")

table(Cluster_group[,2], Cluster_group[,3])
table(Cluster_group[,2], Cluster_group[,4])
table(Cluster_group[,3], Cluster_group[,4])


Cluster_pc_2_dx <- data.frame('Group' = as.factor(Cluster_group$Group), Cluster_pc$x[,1:2]) 
ggplot(data = Cluster_pc_2_dx) + 
  geom_point(aes(x = PC1, y = PC2, col = Group, shape = Group), size = 2) + 
  theme_minimal() 

Cluster_pc_2 <- data.frame('Group' = as.factor(Cluster_group[,2]), Cluster_pc$x[,1:2]) 
ggplot(data = Cluster_pc_2) + 
  geom_point(aes(x = PC1, y = PC2, col = Group, shape = Group), size = 2) + 
  theme_minimal() 

Cluster_pc_3 <- data.frame('Group' = as.factor(Cluster_group[,3]), Cluster_pc$x[,1:2]) 
ggplot(data = Cluster_pc_3) + 
  geom_point(aes(x = PC1, y = PC2, col = Group, shape = Group), size = 2) + 
  theme_minimal() 



Cluster_analysis <- merge(Cluster_group, All_wide, by.x = "V1", by.y = "ID") 
Cluster_analysis[,c("V2")] <- as.factor(as.numeric(Cluster_analysis[,c("V2")]))
Cluster_analysis[,c("V3")] <- as.factor(as.numeric(Cluster_analysis[,c("V3")]))

describe(Cluster_analysis[,c("V2", "V3", "hamd_use", "ymrs_use", "panss_pos_use", "panss_neg_use", "panss_gen_use")])
describeBy(Cluster_analysis[,c("V2","hamd_use", "ymrs_use", "panss_pos_use", "panss_neg_use", "panss_gen_use", "scid_bpi_pf.x")], group = "V2")
describeBy(Cluster_analysis[,c("V3", "hamd_use", "ymrs_use", "panss_pos_use", "panss_neg_use", "panss_gen_use", "scid_bpi_pf.x")], group = "V3")


t.test(hamd_use ~ V2, data=Cluster_analysis, var.equal=TRUE)
t.test(ymrs_use ~ V2, data=Cluster_analysis, var.equal=TRUE)
t.test(panss_pos_use ~ V2, data=Cluster_analysis, var.equal=TRUE)
t.test(panss_neg_use ~ V2, data=Cluster_analysis, var.equal=TRUE)
t.test(panss_gen_use ~ V2, data=Cluster_analysis, var.equal=TRUE)
chisq.test(Cluster_analysis$scid_bpi_pf.x, Cluster_analysis$V2)


anova(aov(hamd_use ~ V3, data=Cluster_analysis))
anova(aov(ymrs_use ~ V3, data=Cluster_analysis))
anova(aov(panss_pos_use ~ V3, data=Cluster_analysis))
anova(aov(panss_neg_use ~ V3, data=Cluster_analysis))
anova(aov(panss_gen_use ~ V3, data=Cluster_analysis))
chisq.test(Cluster_analysis$scid_bpi_pf.x, Cluster_analysis$V3)
chisq.test(Cluster_analysis[(Cluster_analysis$V3 == 1 | Cluster_analysis$V3 == 2), c("scid_bpi_pf.x")], Cluster_analysis[(Cluster_analysis$V3 == 1 | Cluster_analysis$V3 == 2),c("V3")])
chisq.test(Cluster_analysis[(Cluster_analysis$V3 == 1 | Cluster_analysis$V3 == 3), c("scid_bpi_pf.x")], Cluster_analysis[(Cluster_analysis$V3 == 1 | Cluster_analysis$V3 == 3),c("V3")])
chisq.test(Cluster_analysis[(Cluster_analysis$V3 == 2 | Cluster_analysis$V3 == 3), c("scid_bpi_pf.x")], Cluster_analysis[(Cluster_analysis$V3 == 2 | Cluster_analysis$V3 == 3),c("V3")])



cor(Cluster_analysis$scid_bpi_pf.x, Cluster_analysis$panss_gen_use)

PANSS_3Cluster <- aov(panss_gen_use ~ V3, data=Cluster_analysis)
emmeans(PANSS_3Cluster, pairwise ~ "V3", adjust="none")


## Calculation of Intralass correlation to determine the invariance of tasks
# HC
HC_ERN_icc<-lmerTest::lmer(erp_amp ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_ERN[(All_ERN$Group == "HC" & All_ERN$resp_type == "ERN"),])
var_HC_ERN <- as.data.frame(VarCorr(HC_ERN_icc))
# HC ERN Task
var_HC_ERN$vcov[var_HC_ERN$grp == "task_val"] / (var_HC_ERN$vcov[var_HC_ERN$grp == "task_val"] + var_HC_ERN$vcov[var_HC_ERN$grp == "ID"] + var_HC_ERN$vcov[var_HC_ERN$grp == "Residual"])
# HC ERN ID
var_HC_ERN$vcov[var_HC_ERN$grp == "ID"] / (var_HC_ERN$vcov[var_HC_ERN$grp == "task_val"] + var_HC_ERN$vcov[var_HC_ERN$grp == "ID"] + var_HC_ERN$vcov[var_HC_ERN$grp == "Residual"])

HC_CRN_icc<-lmerTest::lmer(erp_amp ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_ERN[(All_ERN$Group == "HC" & All_ERN$resp_type == "CRN"),])
var_HC_CRN <- as.data.frame(VarCorr(HC_CRN_icc))
# HC CRN Task
var_HC_CRN$vcov[var_HC_CRN$grp == "task_val"] / (var_HC_CRN$vcov[var_HC_CRN$grp == "task_val"] + var_HC_CRN$vcov[var_HC_CRN$grp == "ID"] + var_HC_CRN$vcov[var_HC_CRN$grp == "Residual"])
# HC CRN ID
var_HC_CRN$vcov[var_HC_CRN$grp == "ID"] / (var_HC_CRN$vcov[var_HC_CRN$grp == "task_val"] + var_HC_CRN$vcov[var_HC_CRN$grp == "ID"] + var_HC_CRN$vcov[var_HC_CRN$grp == "Residual"])


HC_power_e_icc<-lmerTest::lmer(power ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_power[(All_power$Group == "HC" & All_power$resp_type == "Power_E"),])
var_HC_power_e <- as.data.frame(VarCorr(HC_power_e_icc))
# HC power error Task
var_HC_power_e$vcov[var_HC_power_e$grp == "task_val"] / (var_HC_power_e$vcov[var_HC_power_e$grp == "task_val"] + var_HC_power_e$vcov[var_HC_power_e$grp == "ID"] + var_HC_power_e$vcov[var_HC_power_e$grp == "Residual"])
# HC power error ID
var_HC_power_e$vcov[var_HC_power_e$grp == "ID"] / (var_HC_power_e$vcov[var_HC_power_e$grp == "task_val"] + var_HC_power_e$vcov[var_HC_power_e$grp == "ID"] + var_HC_power_e$vcov[var_HC_power_e$grp == "Residual"])

HC_power_c_icc<-lmerTest::lmer(power ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_power[(All_power$Group == "HC" & All_power$resp_type == "Power_C"),])
var_HC_power_c <- as.data.frame(VarCorr(HC_power_c_icc))
# HC power correct Task
var_HC_power_c$vcov[var_HC_power_c$grp == "task_val"] / (var_HC_power_c$vcov[var_HC_power_c$grp == "task_val"] + var_HC_power_c$vcov[var_HC_power_c$grp == "ID"] + var_HC_power_c$vcov[var_HC_power_c$grp == "Residual"])
# HC power correct ID
var_HC_power_c$vcov[var_HC_power_c$grp == "ID"] / (var_HC_power_c$vcov[var_HC_power_c$grp == "task_val"] + var_HC_power_c$vcov[var_HC_power_c$grp == "ID"] + var_HC_power_c$vcov[var_HC_power_c$grp == "Residual"])


HC_itpc_e_icc<-lmerTest::lmer(itpc ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_itpc[(All_itpc$Group == "HC" & All_itpc$resp_type == "ITPC_E"),])
var_HC_itpc_e <- as.data.frame(VarCorr(HC_itpc_e_icc))
# HC ITPC error Task
var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "task_val"] / (var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "task_val"] + var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "ID"] + var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "Residual"])
# HC ITPC error ID
var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "ID"] / (var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "task_val"] + var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "ID"] + var_HC_itpc_e$vcov[var_HC_itpc_e$grp == "Residual"])

HC_itpc_c_icc<-lmerTest::lmer(itpc ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_itpc[(All_itpc$Group == "HC" & All_itpc$resp_type == "ITPC_C"),])
var_HC_itpc_c <- as.data.frame(VarCorr(HC_itpc_c_icc))
# HC ITPC correct Task
var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "task_val"] / (var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "task_val"] + var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "ID"] + var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "Residual"])
# HC ITPC correct ID
var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "ID"] / (var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "task_val"] + var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "ID"] + var_HC_itpc_c$vcov[var_HC_itpc_c$grp == "Residual"])


# BD
BD_ERN_icc<-lmerTest::lmer(erp_amp ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_ERN[(All_ERN$Group == "BD" & All_ERN$resp_type == "ERN"),])
var_BD_ERN <- as.data.frame(VarCorr(BD_ERN_icc))
# BD ERN Task
var_BD_ERN$vcov[var_BD_ERN$grp == "task_val"] / (var_BD_ERN$vcov[var_BD_ERN$grp == "task_val"] + var_BD_ERN$vcov[var_BD_ERN$grp == "ID"] + var_BD_ERN$vcov[var_BD_ERN$grp == "Residual"])
# BD ERN ID
var_BD_ERN$vcov[var_BD_ERN$grp == "ID"] / (var_BD_ERN$vcov[var_BD_ERN$grp == "task_val"] + var_BD_ERN$vcov[var_BD_ERN$grp == "ID"] + var_BD_ERN$vcov[var_BD_ERN$grp == "Residual"])

BD_CRN_icc<-lmerTest::lmer(erp_amp ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_ERN[(All_ERN$Group == "BD" & All_ERN$resp_type == "CRN"),])
var_BD_CRN <- as.data.frame(VarCorr(BD_CRN_icc))
# BD CRN Task
var_BD_CRN$vcov[var_BD_CRN$grp == "task_val"] / (var_BD_CRN$vcov[var_BD_CRN$grp == "task_val"] + var_BD_CRN$vcov[var_BD_CRN$grp == "ID"] + var_BD_CRN$vcov[var_BD_CRN$grp == "Residual"])
# BD CRN ID
var_BD_CRN$vcov[var_BD_CRN$grp == "ID"] / (var_BD_CRN$vcov[var_BD_CRN$grp == "task_val"] + var_BD_CRN$vcov[var_BD_CRN$grp == "ID"] + var_BD_CRN$vcov[var_BD_CRN$grp == "Residual"])


BD_power_e_icc<-lmerTest::lmer(power ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_power[(All_power$Group == "BD" & All_power$resp_type == "Power_E"),])
var_BD_power_e <- as.data.frame(VarCorr(BD_power_e_icc))
# BD power error Task
var_BD_power_e$vcov[var_BD_power_e$grp == "task_val"] / (var_BD_power_e$vcov[var_BD_power_e$grp == "task_val"] + var_BD_power_e$vcov[var_BD_power_e$grp == "ID"] + var_BD_power_e$vcov[var_BD_power_e$grp == "Residual"])
# BD power error ID
var_BD_power_e$vcov[var_BD_power_e$grp == "ID"] / (var_BD_power_e$vcov[var_BD_power_e$grp == "task_val"] + var_BD_power_e$vcov[var_BD_power_e$grp == "ID"] + var_BD_power_e$vcov[var_BD_power_e$grp == "Residual"])

BD_power_c_icc<-lmerTest::lmer(power ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_power[(All_power$Group == "BD" & All_power$resp_type == "Power_C"),])
var_BD_power_c <- as.data.frame(VarCorr(BD_power_c_icc))
# BD power correct Task
var_BD_power_c$vcov[var_BD_power_c$grp == "task_val"] / (var_BD_power_c$vcov[var_BD_power_c$grp == "task_val"] + var_BD_power_c$vcov[var_BD_power_c$grp == "ID"] + var_BD_power_c$vcov[var_BD_power_c$grp == "Residual"])
# BD power correct ID
var_BD_power_c$vcov[var_BD_power_c$grp == "ID"] / (var_BD_power_c$vcov[var_BD_power_c$grp == "task_val"] + var_BD_power_c$vcov[var_BD_power_c$grp == "ID"] + var_BD_power_c$vcov[var_BD_power_c$grp == "Residual"])


BD_itpc_e_icc<-lmerTest::lmer(itpc ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_itpc[(All_itpc$Group == "BD" & All_itpc$resp_type == "ITPC_E"),])
var_BD_itpc_e <- as.data.frame(VarCorr(BD_itpc_e_icc))
# BD ITPC error Task
var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "task_val"] / (var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "task_val"] + var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "ID"] + var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "Residual"])
# BD ITPC error ID
var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "ID"] / (var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "task_val"] + var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "ID"] + var_BD_itpc_e$vcov[var_BD_itpc_e$grp == "Residual"])

BD_itpc_c_icc<-lmerTest::lmer(itpc ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_itpc[(All_itpc$Group == "BD" & All_itpc$resp_type == "ITPC_C"),])
var_BD_itpc_c <- as.data.frame(VarCorr(BD_itpc_c_icc))
# BD ITPC correct Task
var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "task_val"] / (var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "task_val"] + var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "ID"] + var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "Residual"])
# BD ITPC correct ID
var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "ID"] / (var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "task_val"] + var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "ID"] + var_BD_itpc_c$vcov[var_BD_itpc_c$grp == "Residual"])


# SZ
SZ_ERN_icc<-lmerTest::lmer(erp_amp ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_ERN[(All_ERN$Group == "SZ" & All_ERN$resp_type == "ERN"),])
var_SZ_ERN <- as.data.frame(VarCorr(SZ_ERN_icc))
# SZ ERN Task
var_SZ_ERN$vcov[var_SZ_ERN$grp == "task_val"] / (var_SZ_ERN$vcov[var_SZ_ERN$grp == "task_val"] + var_SZ_ERN$vcov[var_SZ_ERN$grp == "ID"] + var_SZ_ERN$vcov[var_SZ_ERN$grp == "Residual"])
# SZ ERN ID
var_SZ_ERN$vcov[var_SZ_ERN$grp == "ID"] / (var_SZ_ERN$vcov[var_SZ_ERN$grp == "task_val"] + var_SZ_ERN$vcov[var_SZ_ERN$grp == "ID"] + var_SZ_ERN$vcov[var_SZ_ERN$grp == "Residual"])

SZ_CRN_icc<-lmerTest::lmer(erp_amp ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_ERN[(All_ERN$Group == "SZ" & All_ERN$resp_type == "CRN"),])
var_SZ_CRN <- as.data.frame(VarCorr(SZ_CRN_icc))
# SZ CRN Task
var_SZ_CRN$vcov[var_SZ_CRN$grp == "task_val"] / (var_SZ_CRN$vcov[var_SZ_CRN$grp == "task_val"] + var_SZ_CRN$vcov[var_SZ_CRN$grp == "ID"] + var_SZ_CRN$vcov[var_SZ_CRN$grp == "Residual"])
# SZ CRN ID
var_SZ_CRN$vcov[var_SZ_CRN$grp == "ID"] / (var_SZ_CRN$vcov[var_SZ_CRN$grp == "task_val"] + var_SZ_CRN$vcov[var_SZ_CRN$grp == "ID"] + var_SZ_CRN$vcov[var_SZ_CRN$grp == "Residual"])


SZ_power_e_icc<-lmerTest::lmer(power ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_power[(All_power$Group == "SZ" & All_power$resp_type == "Power_E"),])
var_SZ_power_e <- as.data.frame(VarCorr(SZ_power_e_icc))
# SZ power error Task
var_SZ_power_e$vcov[var_SZ_power_e$grp == "task_val"] / (var_SZ_power_e$vcov[var_SZ_power_e$grp == "task_val"] + var_SZ_power_e$vcov[var_SZ_power_e$grp == "ID"] + var_SZ_power_e$vcov[var_SZ_power_e$grp == "Residual"])
# SZ power error ID
var_SZ_power_e$vcov[var_SZ_power_e$grp == "ID"] / (var_SZ_power_e$vcov[var_SZ_power_e$grp == "task_val"] + var_SZ_power_e$vcov[var_SZ_power_e$grp == "ID"] + var_SZ_power_e$vcov[var_SZ_power_e$grp == "Residual"])

SZ_power_c_icc<-lmerTest::lmer(power ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_power[(All_power$Group == "SZ" & All_power$resp_type == "Power_C"),])
var_SZ_power_c <- as.data.frame(VarCorr(SZ_power_c_icc))
# SZ power correct Task
var_SZ_power_c$vcov[var_SZ_power_c$grp == "task_val"] / (var_SZ_power_c$vcov[var_SZ_power_c$grp == "task_val"] + var_SZ_power_c$vcov[var_SZ_power_c$grp == "ID"] + var_SZ_power_c$vcov[var_SZ_power_c$grp == "Residual"])
# SZ power correct ID
var_SZ_power_c$vcov[var_SZ_power_c$grp == "ID"] / (var_SZ_power_c$vcov[var_SZ_power_c$grp == "task_val"] + var_SZ_power_c$vcov[var_SZ_power_c$grp == "ID"] + var_SZ_power_c$vcov[var_SZ_power_c$grp == "Residual"])


SZ_itpc_e_icc<-lmerTest::lmer(itpc ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_itpc[(All_itpc$Group == "SZ" & All_itpc$resp_type == "ITPC_E"),])
var_SZ_itpc_e <- as.data.frame(VarCorr(SZ_itpc_e_icc))
# SZ ITPC error Task
var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "task_val"] / (var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "task_val"] + var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "ID"] + var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "Residual"])
# SZ ITPC error ID
var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "ID"] / (var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "task_val"] + var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "ID"] + var_SZ_itpc_e$vcov[var_SZ_itpc_e$grp == "Residual"])

SZ_itpc_c_icc<-lmerTest::lmer(itpc ~ 1 + (1|task_val) + (1|ID), REML = T, data = All_itpc[(All_itpc$Group == "SZ" & All_itpc$resp_type == "ITPC_C"),])
var_SZ_itpc_c <- as.data.frame(VarCorr(SZ_itpc_c_icc))
# SZ ITPC correct Task
var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "task_val"] / (var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "task_val"] + var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "ID"] + var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "Residual"])
# SZ ITPC correct ID
var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "ID"] / (var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "task_val"] + var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "ID"] + var_SZ_itpc_c$vcov[var_SZ_itpc_c$grp == "Residual"])

