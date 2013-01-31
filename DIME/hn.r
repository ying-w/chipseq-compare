library(DIME)
hn1 = read.table("./data/vs_high1_none1_500_chr1")
hn1[is.na(hn1[,1]),1] = 0
pc = proc.time()
hn1.d1 = DIME(hn1[,1])
print(proc.time()-pc)
save(hn1.d1, file="./data/vs_high1_none1_500_chr1-diff")
#q(save="no")

