#Install Packages ####
BiocManager::install("GenomicAlignments")
install.packages("UniprotR")
install.packages("protti", dependencies = TRUE)
install.packages("r3dmol")
#---
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("seqinr")
install.packages("phangorn")
library(Biostrings) #specify to use "translate" using this package
library(tidyverse)
library(genepop)
library(tidyr)
library(ape)
library(phangorn)
library(dplyr)
library(seqinr) # function "translate overlaps with biostrings

#Set correct working directory####
#getwd()
#setwd("Data/")

#Get DNA sequence####
seq01 <- readDNAStringSet("sequence_01.fasta")

#DNA -> Amino Acid Sequence ####
dna_sequences <- readDNAStringSet("sequence_01.fasta")
#view(dna_sequences)
amino_acid_sequences <- Biostrings::translate(dna_sequences)
#two packages overlapped the same function "translate" so i used "Biostrings::" to specify to use this package in the code above. 
amino_acid_sequences

#Writing the amino acid sequence to a fasta file####
writeXStringSet(dna_sequences, "Seq1AA.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

df <- read.csv("AscNum.csv", header = TRUE)


GetProteinGOInfo(df)


# CLEAN UP #####

# Clear environment

rm(list = ls()) 

# Clear plots
graphics.off()

# Clear packages requires the package pacman to work
p_unload(all)  # Remove all add-ons

# Clear console
cat("\014")  # ctrl+L
  

