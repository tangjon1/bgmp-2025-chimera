library(stringi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(seqinr)
#library(Biostrings)

setwd("/data2/calinsplayground/100824_LL_HK_DcuS/")

#####
dfmapfile = "merged.sorted_noN.fasta"
starcodefile = "merged.barcodes_collapse_d1.tsv"
fileoutfasta = "merged_after_consesnus_R.fasta"
#####

df = read.fasta(dfmapfile, seqtype = "DNA", as.string = TRUE, set.attributes=FALSE)
df <- as.data.frame(do.call(rbind, df))
df <- dplyr::rename(df, "seq"="V1")
df$umi <- rownames(df) #add UMIs as column
rownames(df) <- 1:length(df$umi)

df <- df %>%
  tidyr::separate_wider_delim(umi, delim = ".", names = c("id", "dup"),too_few="align_start") %>%
  dplyr::select(-dup) %>%
  dplyr::rename(umi=id)

umis <- df %>%
  group_by(umi) %>%
  summarise(count=n())

#write.table(umis, "UMIs_for_starcode.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

bcs <- read.table(file = paste(starcodefile,sep=""), sep ='\t', header = FALSE)
colnames(bcs) <- c("bc","numo","collapedBCs")

#x<-3
#load data
bcl <- df
colnames(bcl) <- c("ntseq","bc")
rm(df,umis)

#here we make new rows using the collapsed barcodes
bcs <- bcs %>%
  mutate(totalBCsCollapsed=str_count(collapedBCs, ",")) %>%
  separate_rows(collapedBCs) %>%
  dplyr::rename(tempbc=bc, bc=collapedBCs)

bcl <- bcl %>%
  inner_join(bcs,by="bc") %>%
  dplyr::rename(oldbc=bc, bc=tempbc) %>%
  filter(nchar(as.character(ntseq)) > 3) #sequence must have at least this length

rm(bcs)

#how many times do we see each DNA read for each barcode
#this will be used to find consensus
bclnk_uniq_rd_counts_per_bc <- bcl %>% 
  dplyr::select(-oldbc, -numo, -totalBCsCollapsed) %>% 
  dplyr::group_by_all() %>% 
  dplyr::count() %>% 
  ungroup()

#for each unique BC find the DNA sequence with most reads
bclnk_max_for_bc <- bclnk_uniq_rd_counts_per_bc %>%
  group_by(bc) %>%
  summarise(ms=max(n)) %>%
  dplyr::rename(n=ms) %>%
  ungroup()

#use this to filter with a join
#this only keep the majority read for each BC
bcl_max_rd <- bclnk_uniq_rd_counts_per_bc %>%
  semi_join(bclnk_max_for_bc,by=c("bc","n"))

rm(bclnk_uniq_rd_counts_per_bc, bclnk_max_for_bc)

#now we need to address the issue when we have a tie in number of reads and couldn't call majority

#how many reads does each bc have
bcl_max_rd_num <- bcl_max_rd %>% 
  group_by(bc) %>%
  dplyr::count() %>%
  ungroup()


#bcl_mx_nocollisions <- bcl_max_rd_num %>%
#  filter(n==1) %>%
#  select(-n) %>%
#  semi_join(bcl_max_rd,by="bc")

#only one, no collisions (no ties)
bcl_mx_nocoll <- semi_join(bcl_max_rd, 
                                 bcl_max_rd_num %>%
                                   filter(n==1) %>% dplyr::select(-n),
                                 by="bc")


#more than one, tied, collision
bcl_mx_coll <- semi_join(bcl_max_rd, 
                               bcl_max_rd_num %>%
                                 filter(n>1) %>%
                           dplyr::select(-n),
                               by="bc")

rm(bcl_max_rd_num)

#how many
bcl_mx_coll %>% ungroup() %>% dplyr::select(bc) %>% distinct() %>% nrow()

#save bc_join
write.csv(bcl_mx_nocoll, 
          paste(fileoutfasta,"merged_after_consesnus_R_all_info.csv", sep=""), 
          row.names = FALSE, quote=FALSE)

# #this function is used to write out Fasta files
# writeFasta<-function(data, filename){
#   fastaLines = c()
#   for (rowNum in 1:nrow(data)){
#     fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"bc"], sep = "")))
#     fastaLines = c(fastaLines,as.character(data[rowNum,"ntseq"]))
#   }
#   fileConn<-file(filename)
#   writeLines(fastaLines, fileConn)
#   close(fileConn)
# }
# 
# writeFasta(bcl_mx_nocoll, "merged_after_consesnus_R_all_info.fasta")

# Open a file connection to write to the FASTA file
con <- file(fileoutfasta, "w")

# Loop through each row of the data frame and write to the file
for (i in seq_len(nrow(bcl_mx_nocoll))) {
  writeLines(paste0(">", bcl_mx_nocoll$bc[i]), con)
  writeLines(as.character(bcl_mx_nocoll$ntseq[i]), con)
}

# Close the file connection
close(con)


bcl_mx_nocoll <- bcl_mx_nocoll %>%
  mutate(ntlen=nchar(as.character(ntseq)))

ggplot(bcl_mx_nocoll,aes(x=ntlen))+ #%>% filter(ntlen<170 & ntlen>110 )
  geom_histogram(bins=20)+
  xlab("Length (bp)")+
  xlim(565,585)+
  ylab("Counts")

ggsave("plots/lengths.png",width=9,height=6)

bcl_summ <- bcl_mx_nocoll %>%
  group_by(ntlen) %>%
  summarise(count=n())
