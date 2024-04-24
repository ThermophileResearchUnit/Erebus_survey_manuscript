#  Follow the tutorial https://benjjneb.github.io/dada2/tutorial.html

# If you don't have the dada2 library
# https://benjjneb.github.io/dada2/dada-installation.html


# #Load library and check the version
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("dada2")


#install.packages("devtools")
library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions


library(dada2)
packageVersion("dada2")


# Define the path
path <- "C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/DADA2_Analysis_011221/01_Trimmed/filtered/"


# List files to see their names and check you are in the right directory
list.files(path)

fn <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# 
# for (i in fn){
#   gunzip(i)
#   
# }



#fn <- sort( grep(list.files(path, pattern=".fastq", full.names = TRUE), pattern='.gz', invert=TRUE, value=TRUE))


# Read the names in the fastaq files 

# Check quality score of reads ----
plotQualityProfile(fn[1:2])



# plots don't look so good
# truncLen to 265 nt (in spite of "longer" Reverse primer let's be conservative)
# also check for primers


# Check for priemrs ----

# Load the libraries to manipulate DNA strings
library(ShortRead)
packageVersion("ShortRead")

library(Biostrings)
packageVersion("Biostrings")


# Set up the primer sequences 
# where to check primer sequences - https://www.earthmicrobiome.org/protocols-and-standards/16s/

FWD <- "GTGYCAGCMGCCGCGGTAA"   ## 515F (Parada) 
REV <- "CCGYCAATTYMTTTRAGTTT"  ## 926R (Quince) 


# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Primer orientations - forward / complement / reverse / reverse complement
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fns) {
  nhits <- vcountPattern(primer, sread(readFastq(fns)), fixed = FALSE)
  return(sum(nhits > 0))
}






#SKIP THIS!------------

# See if sequences have primer - in this exampla sample #1, change to check others
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn[[1]]))



#                  Forward Complement Reverse RevComp
# FWD.ForwardReads  106305          0       0       0
# REV.ForwardReads       0          0       0   29547

# We need to remove the primers

# How can I remove primers/adapters/etc from my amplicon reads? - check FAQ
# https://benjjneb.github.io/dada2/faq.html

# If primers are at the start of your reads and are a constant length (it's our case) 
# you can use the trimLeft argument of dada2’s filtering functions to remove the primers


# Filter and trim ----


# Give the sample names

sample.names <- sapply(strsplit(basename(fn), "_"), `[`, 3)


# Create a directory for the filtered files
filt <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

out <- filterAndTrim(fn, filt, truncLen=265, trimLeft=19, trimRight = 20,
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=TRUE)
out


# Inspect the quality and the primers
path <- "C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data/filtered/"
# 




#CONTINUE FROM HERE!

path <- "C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/DADA2_Analysis_011221/01_Trimmed/filtered/"




# Read the names in the cutadapt trimmed fastaq.gz files
fn.filt <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Look at primer detection for the a number of samples
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn.filt[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn.filt[[1]]))

#                  Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# REV.ForwardReads       0          0       0       0






plotQualityProfile(fn.filt[1:2])

# Where are the red lines from the prior plot?
# The primers look a lot better, but the quality scores still look appalling


# Learn the errors ----

# Learn the error Rates

#DOESNT WORK
err <- learnErrors(fn.filt, multithread=F)

# visualize the estimated error rates:
plotErrors(err, nominalQ=TRUE)

filt <- fn.filt

# Sample Inference ----

# dada command to get the sequences
dada <- dada(filt, err, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

# 
# Sample 1 - 47574 reads in 11072 unique sequences.
# Sample 2 - 46708 reads in 13832 unique sequences.
# Sample 3 - 41164 reads in 9437 unique sequences.
# Sample 4 - 46789 reads in 11234 unique sequences.
# Sample 5 - 43988 reads in 15947 unique sequences.
# Sample 6 - 58713 reads in 18495 unique sequences.
# Sample 7 - 54838 reads in 11386 unique sequences.
# Sample 8 - 50050 reads in 16211 unique sequences.
# Sample 9 - 51533 reads in 19101 unique sequences.
# Sample 10 - 62291 reads in 12601 unique sequences.
# Sample 11 - 38215 reads in 11027 unique sequences.
# Sample 12 - 31429 reads in 10845 unique sequences.
# Sample 13 - 58228 reads in 15810 unique sequences.
# Sample 14 - 51731 reads in 17201 unique sequences.
# Sample 15 - 16078 reads in 6350 unique sequences.
# Sample 16 - 62766 reads in 19508 unique sequences.
# Sample 17 - 31081 reads in 12518 unique sequences.
# Sample 18 - 52324 reads in 16539 unique sequences.
# Sample 19 - 39086 reads in 9213 unique sequences.
# Sample 20 - 42571 reads in 10784 unique sequences.
# Sample 21 - 49191 reads in 15120 unique sequences.
# Sample 22 - 51087 reads in 13366 unique sequences.
# Sample 23 - 43496 reads in 7940 unique sequences.
# Sample 24 - 49956 reads in 17460 unique sequences.
# Sample 25 - 49331 reads in 10258 unique sequences.
# Sample 26 - 69866 reads in 18189 unique sequences.
# Sample 27 - 55314 reads in 16738 unique sequences.
# Sample 28 - 45502 reads in 14037 unique sequences.
# Sample 29 - 28305 reads in 9150 unique sequences.
# Sample 30 - 20781 reads in 6739 unique sequences.
# Sample 31 - 26397 reads in 7384 unique sequences.
# Sample 32 - 32884 reads in 9443 unique sequences.
# Sample 33 - 47569 reads in 15480 unique sequences.
# Sample 34 - 53591 reads in 12480 unique sequences.
# Sample 35 - 43045 reads in 11426 unique sequences.
# Sample 36 - 50612 reads in 19788 unique sequences.
# Sample 37 - 47584 reads in 11508 unique sequences.
# Sample 38 - 73198 reads in 18338 unique sequences.
# Sample 39 - 49188 reads in 16304 unique sequences.
# Sample 40 - 56463 reads in 17747 unique sequences.
# Sample 41 - 41632 reads in 10504 unique sequences.
# Sample 42 - 38033 reads in 8000 unique sequences.
# Sample 43 - 49095 reads in 13454 unique sequences.
# Sample 44 - 50081 reads in 15815 unique sequences.
# Sample 45 - 43514 reads in 14282 unique sequences.
# Sample 46 - 38171 reads in 10837 unique sequences.
# Sample 47 - 60699 reads in 18271 unique sequences.
# Sample 48 - 33738 reads in 7551 unique sequences.
# Sample 49 - 2514 reads in 855 unique sequences.
# Sample 50 - 33290 reads in 7526 unique sequences.
# Sample 51 - 43940 reads in 13387 unique sequences.
# Sample 52 - 49054 reads in 15729 unique sequences.
# Sample 53 - 108614 reads in 23926 unique sequences.
# Sample 54 - 59749 reads in 15131 unique sequences.
# Sample 55 - 31941 reads in 7651 unique sequences.
# Sample 56 - 77454 reads in 17847 unique sequences.
# Sample 57 - 50138 reads in 14786 unique sequences.
# Sample 58 - 39645 reads in 8577 unique sequences.
# Sample 59 - 36987 reads in 9062 unique sequences.
# Sample 60 - 67104 reads in 18406 unique sequences




# HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32 - recommend for denoising pyrosequencing data (like IT and 454)
# https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data


# Construct sequence table
seqtab <- makeSequenceTable(dada)

dim(seqtab)
# [1]   60 5236






# The sequence table is a matrix with rows corresponding to (and named by) the samples
# and columns corresponding to (and named by) the sequence variants

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# 246 - sequence length
# 427 

# 215 
# 5236 



# Remove chimeras - Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 1905 bimeras out of 5236 input sequences.

dim(seqtab.nochim)
#[1]   60 3331

sum(seqtab.nochim)/sum(seqtab)
# [1] 0.9709162 - the frequency of non-chimeric sequences

# chimeras account for ca 10% of the abundance





# Track reads through the pipeline - look at the number of reads that made it through each step in the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dada, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "dada", "nonchim")
rownames(track) <- sample.names
track

# input filtered   dada nonchim
# 0101  47574        0  47265   46325
# 0102  46708        0  46118   44097
# 0103  41164        0  40891   40007
# 0104  46789        0  46579   46139
# 0105  43988        0  43203   42002
# 0106  58713        0  57838   53685
# 0107  54838        0  54502   52387
# 0108  50050        0  49342   47824
# 0109  51533        0  50255   47618
# 0110  62291        0  61765   57523
# 0111  38215        0  37840   37795
# 0112  31429        0  30695   29336
# 0113  58228        0  57425   53067
# 0114  51731        0  50769   48065
# 0115  16078        0  15579   15017
# 0116  62766        0  61534   58527
# 0117  31081        0  30197   29532
# 0118  52324        0  51402   47932
# 0119  39086        0  38755   38474
# 0120  42571        0  42205   41588
# 0121  49191        0  48274   45551
# 0122  51087        0  50829   50763
# 0123  43496        0  43339   43303
# 0124  49956        0  49161   49110
# 0125  49331        0  49027   48917
# 0133  69866        0  68952   63965
# 0134  55314        0  54425   51578
# 0135  45502        0  44642   43248
# 0136  28305        0  27990   27161
# 0139  20781        0  20436   20143
# 0140  26397        0  26090   25782
# 0142  32884        0  32273   29714
# 0143  47569        0  46546   43990
# 0144  53591        0  53249   53199
# 0146  43045        0  42804   42756
# 0147  50612        0  49126   45425
# 0148  47584        0  47241   45608
# 0149  73198        0  72558   69888
# 0151  49188        0  48667   46714
# 0152  56463        0  55683   54188
# 0153  41632        0  41442   39918
# 0154  38033        0  37748   37680
# 0155  49095        0  48714   48472
# 0156  50081        0  49402   48836
# 0158  43514        0  42982   42932
# 0160  38171        0  37457   35834
# 0161  60699        0  60179   60101
# 0163  33738        0  33558   33543
# 0166   2514        0   2462    2462
# 0167  33290        0  33089   32987
# 0172  43940        0  43373   43228
# 0173  49054        0  48405   48377
# 0174 108614        0 107818  105524
# 0175  59749        0  59213   57302
# 0179  31941        0  31751   31619
# 0181  77454        0  76867   75143
# 0184  50138        0  49394   48292
# 0301  39645        0  39446   39446
# 0303  36987        0  36578   36476
# 0304  67104        0  66240   64284

# outside of filtering there should be no step in which a majority of reads are lost
# have at least thousands of sequences per-sample, to have good resolution down to 1% frequency - https://github.com/benjjneb/dada2/issues/232
# keep a good proportion of the reads (i.e. >>10%) all the way through the pipeline

# Save the ASV table
#MAKE PATH!

setwd("C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/DADA2_Analysis_011221/")
saveRDS(seqtab.nochim, "seqtab.nochim.erebus")

# Save as csv
#MAKE PATH
write.csv(seqtab.nochim, "asv_erebus.csv")

#install.packages("seqinr", repos="http://R-Forge.R-project.org")

#install.packages(seqinr)
# Create a fasta table
library(seqinr)

name.seq <- paste0("ASV", seq(ncol(seqtab.nochim)))
uniqueSeqs <- as.list(colnames(seqtab.nochim))
write.fasta(uniqueSeqs, name.seq, "C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/outputs/uniqueSeqs_erebus.fasta")






# Assign taxonomy ----

# DADA2 provides a native implementation of the naive Bayesian classifier method for this
# Maintains formatted training fastas - https://benjjneb.github.io/dada2/training.html


# Using the naive bayesian classifier ----

# Silva database
taxa.silva <- assignTaxonomy(seqtab.nochim, "C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/silva_nr99_v138_train_set.fa.gz", multithread=F)


#C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/


# Extension
# make species level assignments based on exact matching between ASVs and sequenced reference strains
# Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign species to 16S gene fragments
# paper by Edgar https://academic.oup.com/bioinformatics/article/34/14/2371/4913809

taxaSp.silva <- addSpecies(taxa.silva, "C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/silva_species_assignment_v138.fa.gz", verbose=TRUE) #, multithread=TRUE

# 54 out of 3331 were assigned to the species level.
# Of which 45 had genera consistent with the input table.

# Save the taxanomy table
setwd("C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data/")
saveRDS(taxaSp.silva, "taxaSp.silva")
write.csv(taxaSp.silva, "taxa_erebus_Silva138.csv")




# Tree ----

# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DECIPHER")
# 
# devtools::install_github("KlausVigo/phangorn")


library(DECIPHER)
packageVersion("DECIPHER") #'2.22.0'
library(phangorn)
packageVersion("phangorn") #'2.8.1'




# Get the seqs from seqtab.nochim and store the them so they can be the tips of the tree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs 

s# Align sequences (using DECIPHER)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA,verbose = FALSE)

# Build a tree (using phangorn)
phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

setwd("C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/DADA2_Analysis_011221/")
saveRDS(fitGTR, "fitGTR")

#THIS HAS RUN NOW!






# FINISH THIS!





# Handoff to phyloseq ----

# Load the metadata and environmental data, remove columns with zero data

library(readr)


#envdata <-  read.csv("C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/sample_and_icpms.csv2", sep = ",", row.names = 1, header = TRUE)
envdata <- read.csv("C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/sample_and_icpms2.csv", sep= ",", row.names = 1, header = TRUE)


# Change the row names of the ASV table to match the environmental table 
rownames(envdata)
rownames(seqtab.nochim) <- c("IonCode106","IonCode110")

#"IonCode106","IonCode110"




#install.packages("ape")
library("ape")
# 
# Warning: cannot remove prior installation of package 'ade4'
# Warning in file.copy(savedcopy, lib, recursive = TRUE) :
#   problem copying C:\Users\hr\Documents\R\win-library\4.1\00LOCK\ade4\libs\x64\ade4.dll to C:\Users\hr\Documents\R\win-library\4.1\ade4\libs\x64\ade4.dll: Permission denied
# Warning: restored 'ade4'
# 


# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")


library(phyloseq)
packageVersion("phyloseq") #'1.38.0'

physeq.erebus <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                          tax_table(taxaSp.silva),
                          sample_data(envdata),
                          phy_tree(fitGTR$tree)) #                    sample_data(envdata),



#RENAME THE SAMPLE_DATA BEFORE YOU CAN PUT IT INTO YOUR PHYLOSEQ OBJECT!
saveRDS(physeq.erebus, "C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/physeq.erebus")





# WHEN YOU OPEN THIS FILE ----
physeq.erebus<- readRDS("C:/Users/trine/OneDrive - Danmarks Tekniske Universitet/DTU2/NZ Project 2022/Data_subset/physeq.erebus")

otu_table(physeq.erebus) # it still has the seqs instead of ASV_1, ASV_2, etc...


new.names <- paste0("ASV_", seq(ntaxa(physeq.erebus))) # define new names ASV_1, ASV_2, ...
seqs <- taxa_names(physeq.erebus) # store sequences
names(seqs) <- new.names # make map from ASV to full sequence
taxa_names(physeq.erebus) <- new.names # rename 

# to convert back to the full sequences use "seqs" e.g. seqs["ASV_"] will return the full sequence that was renamed ASV_1



