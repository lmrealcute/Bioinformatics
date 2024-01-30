pacman::p_load(pacman, tidyverse,msa,tidyr,dplyr,genepop,PopGenome)
#PackageDescriptions####
#pacman - loading and unloading packages
#tidyverse - data import, manipulation, & visualization
#Coin - two sample permutation tests 
#FSA - non-parametric posts test
#lmPerm - flexible permutation tests
#Permuco- post-hoc test
#vcdExtra - contingency tables and creating mosaic plots
#psych - some weird stuff (PCA and Factor Analysis)
#rms - regression modeling; "import' function instead of read_csv
#lme4 - linear mixed effect models
#msa
#tidyr
#dplyr
#genepop
#PopGenome
# Plotting Theme #######
mytheme = theme (plot.title = element_text (face = "bold.italic",
                                            size = 20,
                                            hjust = 0.5),
                 plot.subtitle = element_text (face = "italic",
                                               size = 12,
                                               hjust = 0.5),
                 axis.title = element_text(face = "bold",
                                           size = 18),
                 axis.text = element_text(size = 16),
                 legend.title=element_text(size=18),
                 legend.text=element_text(size=16),
                 panel.background=element_blank(),
                 plot.background=element_blank(),
                 panel.grid.major.y=element_line(color="grey",
                                                 linewidth = 0.5))

#import Data in R####
df = read.csv("Data/bird_morphometrics.csv")
#view(df)

summary(df)

#Edit Variables####
df=df%>%
  mutate(Scientific_Name=as.factor(Scientific_Name),
         English_Name = as.factor(English_Name))

summary(df)
attach(df)

#another way to view files####
#view(t(df))
#--
#loading packages####
#lbrary(tidyr)

par(mfrow = c(1,2))
#run code above and then rerun plot codes and then you can compare side by side
#edits parameters
#1 row
# 2 columns 
#two graphs next to eachother 

#  plot it using the 'plot' function
# we use the dollar sign to access column names in the data
# data of this type is called a 'data frame' (columns and rows of data)
plot(df$Wing_Chord, df$Tail_Length)

# we could plot this with log-transformed data
plot(log(df$Wing_Chord), 
     log(df$Tail_Length))


#---
warning
# CLEAN UP #####

# Clear environment

rm(list = ls()) 

# Clear plots
graphics.off()

# Clear packages requires the package pacman to work
p_unload(all)  # Remove all add-ons

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)