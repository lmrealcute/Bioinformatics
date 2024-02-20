#load Packages####
pacman::p_load(pacman, tidyverse, BiocManager, Biostrings, GenomicAlignments, UniprotR, protti, r3dmol, seqinr,phangorn, genepop, tidyr, ape, dplyr,msa)

#Set correct working directory####
getwd()
setwd("Data/")

#Get DNA sequence####
seq <- readDNAStringSet("sequences.fasta")

#Translating DNA -> Amino Acid Sequence ####
dna_sequences <- readDNAStringSet("sequences.fasta")
#view(dna_sequences)
amino_acid_sequences <- Biostrings::translate(dna_sequences)
#two packages overlapped the same function "translate" so i used "Biostrings::" to specify to use this package in the code above. 
#print the sequences to make sure it worked#
amino_acid_sequences
#AAStringSet object of length 20:
#width seq                                                            names               
#[1]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_1
#[2]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_2
#[3]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_3
#[4]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_4
#[5]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_5
#...   ... ...
#[16]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_16
#[17]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_17
#[18]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_18
#[19]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_19
#[20]   214 NSTPRSREGRSQGWA*KSGQSHLLLTFASD...PENFRVSLWDA*CFLSPSFLWLNSCHRKG Homo_sapiens_20

#Writing the amino acid sequence to a fasta file####
#writeXStringSet(amino_acid_sequences, "SeqAA.fasta", append=FALSE,
               #compress=FALSE, compression_level=NA, format="fasta")

#create a consensus####
alignment_set <- DNAStringSet(seq)
consensus <- consensusString(alignment_set)
print(consensus)

#[1] "AACTCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATAAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAACTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAATTCATGTCATAGGAAGGGG"


#Creating an msa alignment####
#install.packages(msa)
#library(msa)

myFirstAlignment <- msa(seq)
myFirstAlignment

print(myFirstAlignment, show="complete")

#convert alignment to seqinR format ####
UnkHumgene <- msa(seq)
UnkHumgene


UnkHumgene <- msaConvert(UnkHumgene, type="seqinr::alignment")
#compute distance matrix using seqinr 
#assuming you want to compute identity distance

d <- dist.alignment(UnkHumgene)

print(d)

#print phylogenetic tree###
UnkHumgene <- nj(d)
plot(UnkHumgene, main = "Phylogenetic Tree of Unknown Human Gene Sequences")

#2 The samples that are more different are 6, 10, and 4. 6 being the most different than the others. 
#3 the gene seems to be a truncated hemoglobin beta chain, the E val is 2e-09. Accession # is KAI2558340

#4 Tranlating sample 6 into an amino acide sequence and writing it to a fasta file####
#Get DNA sequence
seq6 <- readDNAStringSet("Seq6.txt")

#Translating DNA -> Amino Acid Sequence
dna_sequences <- readDNAStringSet("Seq6.txt")
#view(dna_sequences)
amino_acid_sequences <- Biostrings::translate(dna_sequences)
#two packages overlapped the same function "translate" so i used "Biostrings::" to specify to use this package in the code above. 
#print the sequences to make sure it worked#
amino_acid_sequences

#Writing the amino acid sequence to a fasta file
#writeXStringSet(amino_acid_sequences, "Seq6AA.fasta", append=FALSE,
#compress=FALSE, compression_level=NA, format="fasta")

#5 the gene seems to be a truncated hemoglobin beta chain, the E val is 2e-09. Accession # is KAI2558340

#6 diseases associated####
#the following are a list of associated diseases that I found on OMIM
#thalassemia, erythrocytosis, heinz body anemia, sickle cell disease
#due to the phylogenetic tree showing how different sample 6 is and how common sickle cell anemia is. it is highy likely this individual has sickle cell. 

#3D protein structure is in output folder 



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
