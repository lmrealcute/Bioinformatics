#Install Packages ####
if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomicAlignments")
#install.packages("UniprotR")
#install.packages("protti", dependencies = TRUE)
#install.packages("r3dmol")
#install.packages("seqinr")
#install.packages("phangorn")
#library(Biostrings) #specify to use "translate" using this package
#library(tidyverse)
#library(genepop)
#library(tidyr)
#library(ape)
#library(phangorn)
#library(dplyr)
#library(seqinr) # function "translate overlaps with biostrings

#load Packages####
pacman::p_load(pacman, tidyverse, BiocManager, Biostrings, GenomicAlignments, UniprotR, protti, r3dmol, seqinr,phangorn, genepop, tidyr, ape, dplyr)

#Set correct working directory####
getwd()
setwd("Data/")

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
#4. Read this file into R using the appropriate function ####
accession_numbers<- read.table("AccNumbers.txt")

#5. Sample list of accession numbers ####
accession_numbers <- c("M0QT13", "A0A093QZ57", "A0A7K9CTY9", "A0A7K5KG67", "A0A7L2PPG0")

# Convert the list to a character string
accession_string <- paste(accession_numbers, collapse = ",")

# Print the formatted string
print(accession_string)

#6. Reading accession numbers into GetProteinGOInfo ####
AccessionNumbersGO <- GetProteinGOInfo(accession_numbers)
str(AccessionNumbersGO)

#write into csv file from getproteinGOinfo####
#write.csv(AccessionNumbersGO, "AccessionNumbersGO.csv", row.names = FALSE)

#7. Extract GO terms and their counts from AccessionNumbersGO ####
View(AccessionNumbersGO)

#df <-- read.csv(AccessionNumbersGO.csv) couldnt get to work for some reason
#PlotGoInfo(AccessionNumbersGO) #--> Did not work
go_terms <- unlist(strsplit(AccessionNumbersGO$Gene.Ontology..GO., ";"))
go_terms <- gsub("\\[.*?\\]", "", go_terms)  # Remove GO IDs from GO terms

# Create a data frame with GO terms and their counts
go_counts <- data.frame(GoTerm = go_terms, Count = rep(1, length(go_terms)))

# Summarize the counts for each GO term
go_counts <- aggregate(Count ~ GoTerm, go_counts, sum)

# Plot the GO information
barplot(go_counts$Count, names.arg = go_counts$GoTerm,
        xlab = "GO Terms", ylab = "Count", main = "GO Term Distribution")

#10. Use GetPathology_Biotech() and Get.diseases() to find information on any diseases or pathologies associated with your gene ####
GetPathology_Biotech(accession_numbers)
#NA on all counts
Get.diseases(accession_numbers)
#Error in UseMethod("select") : 
#no applicable method for 'select' applied to an object of class "character"


#11. We are going to access structural information using the protti package ####

viewtibble <- fetch_uniprot(accession_numbers)
View(viewtibble)

#12. Pull any available structural information from the Protein DataBase
fetch_pdb("SCN4A")
#None of the provided IDs could be retrieved!


#13. Get information on any available 3D structures for your gene
fetch_alphafold_prediction(accession_numbers)


WARNING
# CLEAN UP #####

# Clear environment

rm(list = ls()) 

# Clear plots
graphics.off()

# Clear packages requires the package pacman to work
p_unload(all)  # Remove all add-ons

# Clear console
cat("\014")  # ctrl+L
  

