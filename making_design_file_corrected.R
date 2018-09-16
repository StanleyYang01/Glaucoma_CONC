
output_file="./data/design_file.txt"
rna.file = "../conc_libby_new/data/Glaucoma_all_gene_counts.txt"
input_file = "../conc_libby_new/data/design_file_pre.csv"

library(tidyverse)
library(readxl)

data.raw <- read.table(rna.file,  sep="\t", head=T, quote="", check.names=F)
Sample_ID_Full <- names(data.raw)[2:ncol(data.raw)]

design_file_pre <- read.csv(input_file, head=T, check.names=F, stringsAsFactors=F)

design_file_pre$Genotype <- gsub("/", "_", design_file_pre$Genotype)

# correct genotype on mouse ID 42013 to Ddit3
design_file_pre%>%filter(Sample_ID %in% c("42013L", "42013R"))
design_file_pre[5:6, "Genotype"] <- c("Ddit3", "Ddit3") 
design_file_pre%>%filter(Sample_ID %in% c("42013L", "42013R"))

# check if the Sample_ID in design_file_pre matches the Sample_ID from the data.raw table

x <- str_split_fixed(Sample_ID_Full, pattern="_", n=6)
all(x[,1] == design_file_pre[,1])

design_file_pre$Sample_ID_Full <- Sample_ID_Full
names(design_file_pre)

design_file_pre <- design_file_pre %>% mutate(ID_simple = paste(Genotype,Treatment,Sample_ID, sep="."), Group = paste(Genotype,Treatment,sep="."))
## extract mouse ID from Sample_ID column
## next to join the mouse sex info from "Jun.Ddit3.Cohort for RNAseq" using mouse ID as key
design_file_pre<-design_file_pre %>% mutate(mouse_ID = str_sub(Sample_ID, 1,5))
libby_sheet <- read_excel("../conc_libby_new/data/Jun.Ddit3 Cohort for RNA Seq.xlsx", sheet = "Sheet1", col_names = TRUE) %>%
               select(`Animal Number`:Gen, Sex)

data.design <- left_join(design_file_pre, libby_sheet, by = c("mouse_ID"="Animal Number"))

write.table(data.design, output_file, sep="\t", row=F, quote=F)
