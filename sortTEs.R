repmask <- read.table(file = "nanoporeParis.fa.out", header = FALSE, sep ='\t')
colnames(repmask) <- c("SW_score","perc_div.","perc_del.","perc_ins.","query_sequence","position_in_query_begin","position_in_query_end","position_in_query_left","+/-","matching_repeat","repeat_class/family","position_in_repeat_begin","position_in_repeat_end","position_in_repeat_left","ID","*")

paris <- repmask[repmask$`matching_repeat` == "Paris",]
parisNanopore <- paris[5]
print(parisNanopore,row.names = FALSE)
write.table(parisNanopore, 
            file = "paris_nanopore.bed", 
            quote = F, 
            sep = "\t", 
            row.names = F,
            col.names = F)