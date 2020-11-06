#!/usr/bin/env Rscript
# coding: utf-8
# Copyright (C) 2020 Eleanor Fewings
#
# Contact: eleanor.fewings@bioquant.uni-heidelberg.de
#
# ====================
# GNU-GLPv3:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# A full copy of the GNU General Public License can be found on
# http://www.gnu.org/licenses/.
#
# ============================================================
# DESCRIPTION:
# Processes tumour volume data to generate GLM and t.tests for 
# tumour x treatment analysis
# ============================================================
# Usage: (from terminal)
#
# $ ./tumourvol_time.R [options]
#
# Options:
#  -i INPUT, --input=INPUT
#  Path to CSV containing formatted tumour volume data (see example) [required]
# 
#  -c CONTROL, --control=CONTROL
#  Name of control treatment [required]
#  
#  -h, --help
#  Show this help message and exit
#
# Required format for input is (where each day is the change in tumour size):
#
# modelID treatment rep day1 day4 day7 day10 day14 day17 ... dayN
# CN1574   vehicle   1    0  110  150   261   370   527 ...  841
# CN1574   vehicle   2    0   37   41    81   135   320 ...  421
# CN1574   vehicle   3    0   42   62   170   209   358 ...  379
# CN1574   vehicle   4    0   12   74   105   159   212 ...  322
# CN1574   vehicle   5    0   19   65   127   243   248 ...  273
# CN1574   vehicle   6    0   15   79   104   164   221 ...  272
#
# ============================================================

#############
## Startup ##
#############

# Clean up
rm(list=ls())
#ggpubr
#rstatix
# Load libraries
libs <- c("dplyr", "tidyr", "ggplot2", "magrittr", "stringr", "glmmTMB", "effectsize", "optparse", "Rmisc", "rstatix")

for (i in libs) {
  if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) { 
    install.packages(i, repos = "https://ftp.fau.de/cran/")
    if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) {
      stop(paste("Unable to install package: ", i, ". Please install manually and restart.", sep=""))
    }
  }
}

## Get options
option_list <- list(
  make_option(c("--input", "-i"), action="store", default=NULL, type='character',
              help="Path to CSV containing formatted tumour volume data (see example) [required]"),
  make_option(c("--control", "-c"), action="store", default=NULL, type='character',
              help="Name of control treatment [required]")
)

opt <- parse_args(OptionParser(option_list=option_list))

## Check for input and output options
if (is.null(opt$input)) {
  message("ERROR: Input missing, please specify input directory with --input, -i flags.")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

if (is.null(opt$control)) {
  message("ERROR: Control missing, please specify name of control treatment with --control, -c flags.")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}


#Save inputs as normal names
input <- opt$input
control <- opt$control

################
## Read Input ##
################

# import data
data <- read.table(input, header = TRUE, stringsAsFactors = FALSE, sep=",")

#Put data into long format
long <- gather(data, key="time", value="volume", 4:ncol(data))

#Remove , from volumes
long$volume <- str_replace_all(long$volume, ",", "")
long$volume <- as.numeric(long$volume)

#Convert time to factor
long$time <- factor(long$time, levels=unique(long$time))

#Set order of treatments
long$treatment <- factor(long$treatment, levels=unique(long$treatment))

#Add repeat number so each tumour can be tracked individually
long$id <- paste(long$modelID, long$treatment, long$rep, sep="_")

#Remove missing values
long <- long[!is.na(long$volume),]

#Remove + sign in treatment as it is treated as regex
long$treatment <- str_replace_all(long$treatment, "\\+", "_")

#Create cumulative tumour volume data
vol <- long %>% group_by(modelID, treatment, rep) %>% mutate(volsum=cumsum(volume))

############
## Plot 1 ##
############

#Run t.test for each day and treatment vs control

#Summarise cummulative tumour volume over repetitions
sums.1 <- summarySE(vol, measurevar="volsum", groupvars=c("modelID","treatment","time"))

#Set ttest data using change in volume size
ttest.1 <- data.frame(modelID=long$modelID, time=long$time, treatment=long$treatment) %>% unique()

pval <- function(x, y){
  modelID <- x[["modelID"]]
  time <- x[["time"]]
  treatment <- x[["treatment"]]
  # Filter data on above
  treatdat <- y[y$modelID == modelID & y$time == time & y$treatment == treatment,]$volume
  controldat <- y[y$modelID == modelID & y$time == time & y$treatment == control,]$volume
  p.adjust(t.test(controldat, treatdat)$p.val, "fdr", 4)
}

#Calculate p values for t test per model and time point for change in tumour volume
ttest.1$pval <- ttest.1 %>% apply(1, pval, y=long)

#Merge with summary for plotting
sums.1 <- merge(sums.1, ttest.1, by=c("modelID", "time", "treatment"), all=TRUE)

#Replace NaN with NA
sums.1$pval[sums.1$pval == "NaN"] <- NA

#Store plot data
p1.dat <- sums.1

#Format p values for 2 decimal places for plot
sums.1$pval <- formatC(sums.1$pval, format = "e", digits=2)
sums.1$pval[sums.1$pval == " NA" | sums.1$pval == "NA"] <- NA

#Store plot
p1 <- ggplot(sums.1, aes(time, volsum, colour=treatment, group=treatment)) +
  geom_point() +
  geom_line(aes(x=time, y=volsum, colour=treatment)) +
  facet_grid(~modelID) +
  geom_errorbar(aes(ymin=volsum-se, ymax=volsum+se), width=.1) +
  geom_text(aes(label=ifelse(is.na(pval), "", pval)), position="stack", size=4)

#Cleanup
rm(sums.1, ttest.1)

############
## Plot 2 ##
############
#Run t.test for each day and treatment vs control with tumour models combined

#Summarise cummulative tumour volume over repetitions without modelID
sums.2 <- summarySE(vol, measurevar="volsum", groupvars=c("treatment","time"))

#Set ttest data using change in volume size
ttest.2 <- data.frame(time=long$time, treatment=long$treatment) %>% unique()

#Change pval function to remove modelID
pval <- function(x, y){
  time <- x[["time"]]
  treatment <- x[["treatment"]]
  # Filter data on above
  treatdat <- y[y$time == time & y$treatment == treatment,]$volume
  controldat <- y[y$time == time & y$treatment == control,]$volume
  p.adjust(t.test(controldat, treatdat)$p.val, "fdr", 4)
}

#Gather p values per treatment
ttest.2$pval <- ttest.2 %>% apply(1, pval, y=long)

#Merge with summary for plotting
sums.2 <- merge(sums.2, ttest.2, by=c("time", "treatment"), all=TRUE)

#Replace NaN with NA
sums.2$pval[sums.2$pval == "NaN"] <- NA

#Store plot data
p2.dat <- sums.2

#Format p values for 2 decimal places for plot
sums.2$pval <- formatC(sums.2$pval, format = "e", digits=2)
sums.2$pval[sums.2$pval == " NA" | sums.2$pval == "NA"] <- NA

#Store plot
p2 <- ggplot(sums.2, aes(time, volsum, colour=treatment, group=treatment)) +
  geom_point() +
  geom_line(aes(x=time, y=volsum, colour=treatment)) +
  geom_errorbar(aes(ymin=volsum-se, ymax=volsum+se), width=.1) +
  geom_text(aes(label=ifelse(is.na(pval), "", pval)), position="stack", size=4)

#Cleanup
rm(sums.2, ttest.2)

############
## Plot 3 ##
############
#Create data showing deviation of growth rate from emtpy vehicle control

#Summarise change in tumour volume over repetitions
sums.3 <- summarySE(long, measurevar="volume", groupvars=c("modelID","treatment","time"))

#Create day as numeric value
sums.3$day <- str_replace_all(sums.3$time, "day", "") %>% as.numeric()

#Subset vehicle data only
sums.vehicle <- sums.3[sums.3$treatment == control,] %>% select(-2, -4, -(6:8))

colnames(sums.vehicle)[3] <- control

#Merge on vehicle as seperate column per model and time
sums.3 <- merge(sums.3, sums.vehicle, by=c("modelID", "time", "day"), all=TRUE)

#zero centred normalisation on vehicle
sums.3$norm <- sums.3$vehicle - sums.3$volume

#Relevel to set vehicle as reference
sums.3$treatment <- factor(sums.3$treatment)
sums.3$treatment <- relevel(sums.3$treatment, ref = "vehicle")

#Function to create model per tumour model
perTglm <- function(x, y) {
  # Filter data on above
  model.dat <- y[y$modelID == x,]
  #Create model
  glm.dat <- glmmTMB(norm ~ day*treatment + 0, data=model.dat)
  #Create dataframe of model summary
  glm.df <- summary(glm.dat)$coefficients$cond %>% as.data.frame() %>% filter(!grepl("day", row.names(.))) %>% select(-Estimate)
  #Calculate effectsizes
  glm.effect <- effectsize(glm.dat) %>% as.data.frame() %>% filter(!grepl("day", Parameter))
  #Label modelID and treatment
  glm.df$modelID <- x
  #Rename columns of glm results
  colnames(glm.df) <- c("Std.err", "Zvalue", "pval", "modelID")
  #Bind results on
  all.glm <- merge(glm.df, glm.effect, by.x="row.names", by.y="Parameter")
  #Change column 1 name
  colnames(all.glm)[1] <- "treatment" 
  all.glm$treatment <- str_replace_all(all.glm$treatment, "treatment", "")
  #Return
  all.glm
}

all.glm <- sums.3$modelID %>% unique() %>% lapply(perTglm, y=sums.3) %>% do.call("rbind", .)

#Set time for pvalue label
all.glm$time <- "day21"

#Merge results onto summary to plot
sums.3 <- merge(sums.3, all.glm, by=c("modelID", "treatment", "time"), all=TRUE)

#Store plot data
p3.dat <- sums.3

#Format p values for 2 decimal places for plot
sums.3$pval <- formatC(sums.3$pval, format = "e", 2)
sums.3$pval[sums.3$pval == " NA" | sums.3$pval == "NA"] <- NA

#Plot data with glm generated p values
p3 <- ggplot(sums.3, aes(time, volume, colour=treatment, group=treatment)) +
  geom_point() +
  geom_line(aes(x=time, y=volume, colour=treatment)) +
  facet_grid(~modelID) +
  geom_errorbar(aes(ymin=volume-se, ymax=volume+se), width=.1) +
  geom_text(aes(label=ifelse(is.na(pval), "", pval)), position="stack", size=4)

#Cleanup
rm(all.glm, sums.vehicle, sums.3)

############
## Plot 4 ##
############
#Create data showing deviation of growth rate from emtpy vehicle control, combining models - using vehicle as reference

#Summarise change in tumour volume over repetitions excluding modelID
sums.4 <- summarySE(long, measurevar="volume", groupvars=c("treatment", "time"))

#Create day as numeric value
sums.4$day <- str_replace_all(sums.4$time, "day", "") %>% as.numeric()

#Subset vehicle data only
sums.vehicle <- sums.4[sums.4$treatment == "vehicle",] %>% select(-1, -3, -(5:7))

colnames(sums.vehicle)[2] <- control

#Merge on vehicle as seperate column per time
sums.4 <- merge(sums.4, sums.vehicle, by=c("day", "time"), all=TRUE)

#zero centred normalisation on vehicle
sums.4$norm <- sums.4$vehicle - sums.4$volume

#Relevel to set vehicle as reference
sums.4$treatment <- factor(sums.4$treatment)
sums.4$treatment <- relevel(sums.4$treatment, ref = "vehicle")

#Create model comparing to volume
mod.all <- glmmTMB(norm ~ day*treatment + 0, data=sums.4)

#Create dataframe of model summary
all.glm <- summary(mod.all)$coefficients$cond %>% as.data.frame() %>% filter(!grepl("day", row.names(.))) %>% select(-Estimate)

#Calculate effectsizes
glm.effect <- effectsize(mod.all) %>% as.data.frame() %>% filter(!grepl("day", Parameter))

#Rename columns of glm results
colnames(all.glm) <- c("Std.err", "Zvalue", "pval")

#Merge on effect size
all.glm <- merge(all.glm, glm.effect, by.x="row.names", by.y="Parameter")

#Change column 1 name
colnames(all.glm)[1] <- "treatment" 
all.glm$treatment <- str_replace_all(all.glm$treatment, "treatment", "")

#Set time for merging
all.glm$time <- "day21"

#Merge results onto summary to plot
sums.4 <- merge(sums.4, all.glm, by=c("treatment", "time"), all=TRUE)

#Store plot data
p4.dat <- sums.4

#Format p values for 2 decimal places for plot
sums.4$pval <- formatC(sums.4$pval, format = "e", 2)
sums.4$pval[sums.4$pval == " NA" | sums.4$pval == "NA"] <- NA

#Plot data with glm generated p values
p4 <- ggplot(sums.4, aes(time, volume, colour=treatment, group=treatment)) +
  geom_point() +
  geom_line(aes(x=time, y=volume, colour=treatment)) +
  geom_errorbar(aes(ymin=volume-se, ymax=volume+se), width=.1) +
  geom_text(aes(label=ifelse(is.na(pval), "", pval), position="stack", size=4))

#Cleanup
rm(all.glm, sums.vehicle, sums.4)

###################
## Write outputs ##
###################

#Rename columns to make them more intuitive
colnames(p1.dat) <- c("modelID", "time", "treatment", "repeats", "cumulative.volume", "standard.dev", "standard.err", "confidence.int", "adj.pval")
colnames(p2.dat) <- c("time", "treatment", "repeats", "cumulative.volume", "standard.dev", "standard.err", "confidence.int", "adj.pval")
colnames(p3.dat) <- c("modelID", "treatment", "time", "day", "repeats", "cumulative.volume", "standard.dev", "standard.err", "confidence.int", "control.volume", "control.normalised.volume", "glm.standard.err", "glm.Z.Value", "glm.P.value", "glm.standard.coefficient", "glm.CI", "glm.CI.low", "glm.CI.high")
colnames(p4.dat) <- c("treatment", "time", "day", "repeats", "cumulative.volume", "standard.dev", "standard.err", "confidence.int", "control.volume", "control.normalised.volume", "glm.standard.err", "glm.Z.Value", "glm.P.value", "glm.standard.coefficient", "glm.CI", "glm.CI.low", "glm.CI.high")

#Set output directory
out <- dirname(input)

#Write data files
write.csv(p1.dat, paste(out, "P1.tumourvolume.ttestxmodel.dat.csv", sep="/"), row.names = FALSE)
write.csv(p2.dat, paste(out, "P2.tumourvolume.ttest.dat.csv", sep="/"), row.names = FALSE)
write.csv(p3.dat, paste(out, "P3.tumourvolume.glmxmodel.dat.csv", sep="/"), row.names = FALSE)
write.csv(p4.dat, paste(out, "P4.tumourvolume.glm.dat.csv", sep="/"), row.names = FALSE)

#Write plots
ggsave("P1.tumourvolume.ttestxmodel.pdf", p1, width = 16, height=9)
ggsave("P2.tumourvolume.ttest.pdf", p2, width = 16, height=9)
ggsave("P3.tumourvolume.glmxmodel.pdf", p3, width = 16, height=9)
ggsave("P4.tumourvolume.glm.pdf", p4, width = 16, height=9)

