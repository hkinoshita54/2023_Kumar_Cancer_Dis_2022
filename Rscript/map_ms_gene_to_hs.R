library(biomaRt)
library(org.Mm.eg.db)   # for mouse cells
library(tidyverse)

# create a table to map mouse gene names to human symbol
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() 
bm[bm == ""] <- NA
bm <- bm %>% na.omit()

# join to the result table
MUC6KO_UP <- read.delim("gene_set/MUC6KO_UP_500.txt", header = F, sep = " ")
MUC6KO_UP <- inner_join(MUC6KO_UP, bm, by = c("V1"="external_gene_name"))
write.table(MUC6KO_UP,"gene_set//MUC6KO_UP_hs.txt", sep ="\t", col.names = T,row.names = F)

MUC6KO_DN <- read.delim("gene_set/MUC6KO_DN.txt", header = F, sep = " ")
MUC6KO_DN <- inner_join(MUC6KO_DN, bm, by = c("V1"="external_gene_name"))
write.table(MUC6KO_DN,"gene_set//MUC6KO_DN_hs.txt", sep ="\t", col.names = T,row.names = F)

# join to the result table for ARID1A
ARID1AKO_UP <- read.delim("gene_set/ARID1AKO_UP.txt", header = F, sep = " ")
ARID1AKO_UP <- inner_join(ARID1AKO_UP, bm, by = c("V1"="external_gene_name"))
write.table(ARID1AKO_UP,"gene_set//ARID1AKO_UP_hs.txt", sep ="\t", col.names = T,row.names = F)

ARID1AKO_DN <- read.delim("gene_set/ARID1AKO_DN.txt", header = F, sep = " ")
ARID1AKO_DN <- inner_join(ARID1AKO_DN, bm, by = c("V1"="external_gene_name"))
write.table(ARID1AKO_DN,"gene_set//ARID1AKO_DN_hs.txt", sep ="\t", col.names = T,row.names = F)