library(dplyr)
library('tidyr')

#repeatmasker file
#1: Smith-Waterman score of the match, usually complexity adjusted
  #The SW scores are not always directly comparable. Sometimes
  #the complexity adjustment has been turned off, and a variety of
  #scoring-matrices are used.
#2: % substitutions in matching region compared to the consensus
#3: % of bases opposite a gap in the query sequence (deleted bp)
#4: % of bases opposite a gap in the repeat consensus (inserted bp)
#5: name of query sequence
#6: starting position of match in query sequence
#7: ending position of match in query sequence
#8: no. of bases in query sequence past the ending position of match
#9: match is with the Complement of the consensus sequence in the database
#10: name of the matching interspersed repeat
#11: the class of the repeat, in this case a DNA transposon 
  #fossil of the MER2 group (see below for list and references)
#12: no. of bases in (complement of) the repeat consensus sequence 
  #prior to beginning of the match (so 0 means that the match extended 
  #all the way to the end of the repeat consensus sequence)
#13: starting position of match in database sequence (using top-strand numbering)
#14: ending position of match in database sequence
repmask <- read.table(file = "nanoporeParis.fa.out", header = FALSE, sep ='\t')
colnames(repmask) <- c("SW_score","perc_div.","perc_del.","perc_ins.","query_sequence","position_in_query_begin","position_in_query_end","position_in_query_left","+/-","matching_repeat","repeat_class/family","position_in_repeat_begin","position_in_repeat_end","position_in_repeat_left","ID","*")
paris <- repmask[repmask$`repeat_class/family` == "Unspecified",]

bed160 <- read.table(file="nanoParis160.bed", sep = "\t")
colnames(bed160) <- c("Chr/Scaffold","Start","Stop","Nanopore","Q-Map","+/-")

bed9 <- read.table(file="nanoporeParis9.bed", sep = "\t")
colnames(bed9) <- c("Chr/Scaffold","Start","Stop","Nanopore","Q-Map","+/-")


#merge the repeatmasker output with the bwa and final bedfile
merged160 <- merge(bed160, paris,  by.y= "query_sequence", by.x = "Nanopore")
#sort the merged160 data frame by the chromosome then the starting position
merged160 <- merged160[order(merged160$`Chr/Scaffold`,merged160$Start),]
#reorder to be in bedfile format chromosome, start, stop, name
merged160 <- merged160[c(2,3,4,1,5,15)]
#unite the TE to the naopore read ID
merged160 <- merged160 %>% unite("Nanopore_TE",c(4,6),remove = T)
#write to a file to save new bed file
write.table(merged160, 
            file = "merged160.bed", 
            quote = F, 
            sep = "\t", 
            row.names = F,
            col.names = F)


#merge the repeatmasker output with the bwa and final bedfile
merged9 <- merge(bed9, paris,  by.y= "query_sequence", by.x = "Nanopore")
#sort the merged160 data frame by the chromosome then the starting position
merged9 <- merged9[order(merged9$`Chr/Scaffold`,merged9$Start),]
#reorder to be in bedfile format chromosome, start, stop, name
merged9 <- merged9[c(2,3,4,1,5,15)]
merged9 <- merged9 %>% unite("Nanopore_TE",c(4,6),remove = T)
write.table(merged9, 
            file = "merged9.bed", 
            quote = F, 
            sep = "\t", 
            row.names = F,
            col.names = F)

mergeBedRepmask <- function(bedFile,repmaskFile){
  merged <- merge(bedFile, repmaskFile,  by.y= "query_sequence", by.x = "Nanopore")
  merged <- merged[order(merged$`Chr/Scaffold`,merged$Start),]
  merged <- merged[c(2,3,4,1,5,15)]
  merged <- merged %>% unite("Nanopore_TE",c(4,6),remove = T)
  write.table(merged, 
              file = "merged.bed", 
              quote = F, 
              sep = "\t", 
              row.names = F,
              col.names = F)
  print(merged)
}

mergeBedRepmask(bed160,paris)
