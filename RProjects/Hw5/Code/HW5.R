if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("seqinr")
install.packages("phangorn")
library(Biostrings)
library(tidyverse)
library(genepop)
library(tidyr)
library(ape)
library(phangorn)
library(dplyr)
library(seqinr)

setwd("Data/")
#set working directory for the HW####
getwd()
#Check Wd in good####
seq01 <- readDNAStringSet("sequence_01.fasta")
seq02 <- readDNAStringSet("sequence_02.fasta")
seq03 <- readDNAStringSet("sequence_03.fasta")
seq04 <- readDNAStringSet("sequence_04.fasta")
seq05 <- readDNAStringSet("sequence_05.fasta")
combseq <- readDNAStringSet("CombinedFrogSequence.txt")

Combinedseq <- c(seq01, seq02, seq03, seq04, seq05)
Combinedseq

#create a consensus####
alignment_set <- DNAStringSet(Combinedseq)
consensus <- consensusString(alignment_set)
print(consensus)


#Creating an msa alignment####
install.packages(msa)
library(msa)

myFirstAlignment <- msa(Combinedseq)
myFirstAlignment

print(myFirstAlignment, show="complete")

#Function to count gaps in consensus####
if (!require(stringr))
install.packages("stringr") 
library(stringr)

str_count(consensus, "-")
# 0 gaps

#calculate and print the width of alignments####
alignment_length <-width(consensus)
print(alignment_length)
#alignment length is 566 bp 

#function to count GC Content in a concensus####
library(dplyr)
#--

Table <- myFirstAlignment  %>% 
  paste(collapse="") %>% 
  strsplit(split="") %>% unlist %>% 
  `[`(!. %in% c("", " ", ".", ",")) %>% 
  table

A <-  Table["A"]
T <-  Table["T"]
G <-  Table["G"]
C <-  Table["C"]

Total <- A + T + G + C

GC <- (G + C) / Total
GC
#GC Content for consensus is 0.3856631 or 38% 

#convert alignment to seqinR format ####
FrogSCN4Agene <- msa(Combinedseq)
FrogSCN4Agene


FrogSCN4Agene <- msaConvert(FrogSCN4Agene, type="seqinr::alignment")
#compute distance matric using seqinr 
#assuming you want to compute identity distance

d <- dist.alignment(FrogSCN4Agene)

print(d)

#print phylogenetic tree###
FrogSCN4Agene <- nj(d)
plot(FrogSCN4Agene, main = "Phylogenetic Tree of Frog SCN4A Gene Sequences")
#

#Translate DNA sequence into AA sequence ####
dna_sequences <-readDNAStringSet("Data/sequence_01.fasta")
amino_acid_sequences <- translate(dna_sequences)
amino_acid_sequences

#convert alignment to phangorn####
Alignment_phyDat <- msaConvert(myFirstAlignment, type = "phangorn::phyDat")
Alignment_phyDat
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")

