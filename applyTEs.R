library(dplyr)

repmask <- read.table("nanoporeParis.fa.out", sep ='\t')
print(repmask[c(5,10)])

bed160 <- read.table(file="nanoParis160.bed")
