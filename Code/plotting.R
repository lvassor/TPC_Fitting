#!/usr/bin/env Rscript

# Author: Luke Vassor login: ljv3618@ic.ac.uk
# Script: plotting.R
# Desc: R script which plots Thermal Performance Curves using fitting output
# Date: Mar 2019

# Clear workspace
rm(list=ls())
graphics.off()

# Import and Libraries
require(ggplot2, warn.conflicts=FALSE)
library(extrafont)
library(xtable)
library(tidyr, warn.conflicts=FALSE)
loadfonts(quiet = TRUE)
# pdfFonts()

###############################
#####      WRANGLING     ######
###############################
# Read in data
print("Importing model fitting results...")

ODir <- "../Results/"
dir.create(paste(ODir,'Figs', sep = ""), showWarnings = FALSE)
dir.create(paste(ODir,'Figs/Model_Fits', sep = ""), showWarnings = FALSE)
dir.create(paste(ODir,'Figs/Stats', sep = ""), showWarnings = FALSE)

wrangled_data <- as.data.frame(read.csv("../Results/wrangled_data.csv"))
model_fit_results_full <- as.data.frame(read.csv("../Results/model_fit_results.csv"), header = TRUE)
# model_fit_results_full <- na.omit(model_fit_results_full) #[complete.cases(model_fit_results_full), ]

IDs <- unique(wrangled_data$FinalID)
IDs <- as.character(IDs)
names_to_keep <- c("FinalID",
                            "sch_final_B_0", "sch_final_E", "sch_final_E_h", "sch_final_T_h",
                            "bri_final_B_0", "bri_final_T_0", "bri_final_T_m", 
                            "cub_final_B_0", "cub_final_B_1", "cub_final_B_2", "cub_final_B_3",
                            "sch_R_Squared", "sch_AICc",
                            "bri_R_Squared", "bri_AICc",
                            "cub_R_Squared", "cub_AICc",
                            "sch_delta_AICc", "bri_delta_AICc", "cub_delta_AICc",
                            "sch_Wi_AICc", "bri_Wi_AICc", "cub_Wi_AICc", "Points_Before_Peak", "Points_After_Peak")
model_fit_results <- subset(model_fit_results_full, select = names_to_keep)
###############################
#####     STATISTICS     ######
###############################

# Obtaining results - means
results <- matrix(1:9, nrow = 3)
rownames(results) <- c("SSH", "Briere", "Cubic")
colnames(results) <- c("Successful Fits (%)", "Mean R^2", "Mean AICc")
# Calculate % of successful fits (unsuccessful fit = NA)
results[1,1] = 100*(length(model_fit_results[,"sch_final_B_0"]) - sum(is.na(model_fit_results[,"sch_final_B_0"])))/length(model_fit_results[,"sch_final_B_0"])
results[2,1] = 100*(length(model_fit_results[,"bri_final_B_0"]) - sum(is.na(model_fit_results[,"bri_final_B_0"])))/length(model_fit_results[,"bri_final_B_0"])
results[3,1] = 100*(length(model_fit_results[,"cub_final_B_0"]) - sum(is.na(model_fit_results[,"cub_final_B_0"])))/length(model_fit_results[,"cub_final_B_0"])

# Calculate model mean R2 values over all fits
# subset to omit rows with negative R2 (these would effect successful fits above)
model_fit_results <- subset(model_fit_results, sch_R_Squared > 0 & bri_R_Squared > 0 & cub_R_Squared > 0)
results[1,2] = mean(c(model_fit_results[,"sch_R_Squared"]))
results[2,2] = mean(c(model_fit_results[,"bri_R_Squared"]), na.rm = T) # only briere has NA values, adapt as you wish
results[3,2] = mean(c(model_fit_results[,"cub_R_Squared"]))

# Calculate Mean AICc
results[1,3] = mean(c(model_fit_results[,"sch_AICc"]))
results[2,3] = mean(c(model_fit_results[,"bri_AICc"]), na.rm = T)
results[3,3] = mean(c(model_fit_results[,"cub_AICc"]))
# xtable(results)

# Obtaining results - Akaike weights
model_fit_results_all <- na.omit(model_fit_results)
delta_AICc_results <- matrix(1:18, nrow = 3)
rownames(delta_AICc_results) <- c("SSH", "Briere", "Cubic")
colnames(delta_AICc_results) <- c("delta <= 2", "2 < delta <= 4", "4 < delta <= 7", "7 < delta <= 10", "delta > 10", "Total")
schoolfield_delta_AICc <- c(model_fit_results_all[,"sch_delta_AICc"]) # make vectors from columns for analysis
briere_delta_AICc <- c(model_fit_results_all[,"bri_delta_AICc"])
cubic_delta_AICc <- c(model_fit_results_all[,"cub_delta_AICc"])
# Schoolfield
delta_AICc_results[1,1] = sum(schoolfield_delta_AICc <= 2)
delta_AICc_results[1,2] = sum(schoolfield_delta_AICc > 2 & schoolfield_delta_AICc <= 4)
delta_AICc_results[1,3] = sum(schoolfield_delta_AICc > 4 & schoolfield_delta_AICc <= 7)
delta_AICc_results[1,4] = sum(schoolfield_delta_AICc > 7 & schoolfield_delta_AICc <= 10)
delta_AICc_results[1,5] = sum(schoolfield_delta_AICc > 10)
delta_AICc_results[1,6] = length(schoolfield_delta_AICc)
# Briere
delta_AICc_results[2,1] = sum(briere_delta_AICc <= 2)
delta_AICc_results[2,2] = sum(briere_delta_AICc > 2 & briere_delta_AICc <= 4)
delta_AICc_results[2,3] = sum(briere_delta_AICc > 4 & briere_delta_AICc <= 7)
delta_AICc_results[2,4] = sum(briere_delta_AICc > 7 & briere_delta_AICc <= 10)
delta_AICc_results[2,5] = sum(briere_delta_AICc > 10)
delta_AICc_results[2,6] = length(briere_delta_AICc)

# Cubic
delta_AICc_results[3,1] = sum(cubic_delta_AICc <= 2)
delta_AICc_results[3,2] = sum(cubic_delta_AICc > 2 & cubic_delta_AICc <= 4)
delta_AICc_results[3,3] = sum(cubic_delta_AICc > 4 & cubic_delta_AICc <= 7)
delta_AICc_results[3,4] = sum(cubic_delta_AICc > 7 & cubic_delta_AICc <= 10)
delta_AICc_results[3,5] = sum(cubic_delta_AICc > 10)
delta_AICc_results[3,6] = length(cubic_delta_AICc)
# xtable(delta_AICc_results)

# Calculate akaike weights
sch_weights <- c(model_fit_results_all[,"sch_Wi_AICc"])
bri_weights <- c(model_fit_results_all[,"bri_Wi_AICc"])
cub_weights <- c(model_fit_results_all[,"cub_Wi_AICc"])
data <- data.frame(sch_weights, bri_weights, cub_weights)

# Calculate Akaike weights 
###############################
#####      PLOTTING       #####
###############################

# Original data Akaike weights Violin plot
data <- tidyr::gather(data, "Model","Akaike_Weight", 1:3)
weights <- qplot(x=Model, y = Akaike_Weight, data = data, geom = c("violin"), fill = Model) + theme_bw() + theme(text = element_text(family="CM Roman", size = 10))
weights <- weights + scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00"))
weights <- weights + ylab(expression(paste("Akaike Weight, W"[i],"AIC"[c])))
ggsave(filename = paste('../Results/Figs/Stats/', "Akaike_Weights", ".pdf", sep = ""), plot = weights, height = 4, width = 4.2, device = "pdf")
embed_fonts(paste('../Results/Figs/Stats/', "Akaike_Weights", ".pdf", sep = ""), outfile="../Results/Figs/Stats/Akaike_Weights.pdf")


# Define model functions

cubic <- function(B_0, B_1, B_2, B_3, T){
    return(B_0 + B_1*T + B_2*T^2 + B_3*T^3)
}

briere <- function(B_0, T_0, T_m, T){
    return(B_0*T*(T-T_0)*((T_m-T)^0.5))
}

schoolfield <- function(B_0, E, k, T_h, E_h, T){
    return((B_0*exp((-E/k)*((1/T)-(1/283.15))))/(1+exp((E_h/k)*(1/T_h - 1/T))))
}
# Go back to results with negative R2 to produce all plots
model_fit_results <- subset(model_fit_results_full, select = names_to_keep)
start <- proc.time()[3]
print("Beginning TPC plotting...")
for (i in 1:length(IDs)){   # loop to run analysis for each thermal response separately
    # print(i)
    tmp_mdl_values <- subset(model_fit_results, FinalID == IDs[i])
    tmp_trait_temp_values <- subset(wrangled_data, FinalID == IDs[i])
    DataToPlot <- data.frame(Temperature = tmp_trait_temp_values[, "ConTemp_K"], OriginalTraitValue = tmp_trait_temp_values[, "OriginalTraitValue"])
    Trait_Name <- tmp_trait_temp_values[1, "StandardisedTraitName"]
    ## Generate temp sequence for model
    temps <- seq(min(floor(tmp_trait_temp_values$ConTemp_K)), ceiling(max(tmp_trait_temp_values$ConTemp_K)), length = 200)
    ## Schoolfield model

    schoolfield_model <- schoolfield(B_0 = tmp_mdl_values$sch_final_B_0, E = tmp_mdl_values$sch_final_E, E_h = tmp_mdl_values$sch_final_E_h, T_h = tmp_mdl_values$sch_final_T_h, k = 8.617e-05, T = temps)
    ## Briere model
    briere_model <- briere(B_0 = tmp_mdl_values$bri_final_B_0, T_0 = tmp_mdl_values$bri_final_T_0, T_m = tmp_mdl_values$bri_final_T_m, T = temps)
    ## Cubic model
    cubic_model <- cubic(B_0 = tmp_mdl_values$cub_final_B_0, B_1 = tmp_mdl_values$cub_final_B_1, B_2 = tmp_mdl_values$cub_final_B_2, B_3 = tmp_mdl_values$cub_final_B_3, T = temps)

    ## Create plot
    cols <- c("BRI"="#999999", "SSH"="#E69F00", "CUB"="#56B4E9")
    p <- ggplot() + #Plot just the raw data
    geom_point(data = DataToPlot, aes(x = Temperature - 273.15, y = OriginalTraitValue), size = I(2), colour = "black", alpha = 0.9) + 
    xlab(expression(paste("Temperature (",degree,C,")"))) +
    ylab(paste(Trait_Name)) +
    theme_bw() + theme(text = element_text(family="CM Roman", size = 10))

    ## Create Lines
    if (any(is.na(schoolfield_model)) == FALSE){ # mask for each model so don't have to eliminate entire row = more curves
        SchToPlot <- data.frame(Temperature = temps - 273.15, OriginalTraitValue = schoolfield_model)
        p <- p + geom_line(data = SchToPlot, aes(x = Temperature, y = OriginalTraitValue, colour = "SSH"))
    }
    if (any(is.na(briere_model)) == FALSE){
        BriToPlot <- data.frame(Temperature = temps - 273.15, OriginalTraitValue = briere_model)
        p <- p + geom_line(data = BriToPlot, aes(x = Temperature, y = OriginalTraitValue, colour = "BRI"))
    }
    if (any(is.na(cubic_model)) == FALSE){
        CubToPlot <- data.frame(Temperature = temps - 273.15, OriginalTraitValue = cubic_model)
        p <- p + geom_line(data = CubToPlot, aes(x = Temperature, y = OriginalTraitValue, colour = "CUB"))
    }
    ## Add lines
    # p <- p + scale_color_discrete(name = "Model", labels = c("SSH", "Briere", "Cubic Poly"))
    
    p <- p + scale_colour_manual(name="Model", values = cols)
    ggsave(filename = paste('../Results/Figs/Model_Fits/', IDs[i], ".pdf", sep = ""), plot = p, height = 4, width = 4.2, device = "pdf")
    embed_fonts(paste('../Results/Figs/Model_Fits/', IDs[i], ".pdf", sep = ""), outfile=paste('../Results/Figs/Model_Fits/', IDs[i], ".pdf", sep = ""))
}
print("Plots complete! See ../Results/Figs/Model_Fits directory")
# sprintf("Completed in: %f seconds", proc.time()[3] - start[3])