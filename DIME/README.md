# DIME

Sorry this directory is currently a mess.

There is some code here to preprocess data in such a way that [DIME](http://www.stat.osu.edu/~statgen/SOFTWARE/DIME/) can be used.

The author (C. Taslim) was kind enough to send me some matlab code to do the loess normalization (used for input into DIME).
The matlab code can is in the following 3 files:

    normalizeMeanVarStep1.m
    normalizeMeanVarStep2.m
    runNormzMeanVarStep1_2.m
    
I have rewritten parts of what the three matlab scripts will do and put them into an R file: `loess.R`

If memory serves me correctly, I then made non-overlapping bins over an entire chromosome and counted reads in each bin using something like:

```bash
perl covBEDforDIME.pl chromsizes19M.tab > 500bp_window.bed &
#each takes about 6min, do not run in parallel (disk bound, splits into multiple processes)
intersectBed -v -abam ../experimentA1.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprA1.cov 
intersectBed -v -abam ../experimentA2.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprA2.cov 
intersectBed -v -abam ../experimentB1.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprB1.cov
intersectBed -v -abam ../experimentB2.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprB2.cov
```
Where `chromsizes19M.tab` is a two colum file of chromosome names in the first column and chromosome sizes in the second column and 
`artifact2.bed` is merged file containing DAC Blacklisted Regions and Duke Excluded Regions from ENCODE project. 

I then split up the file the coverage files by chromosome using `splitcov.pl` and tried run `loess.R` and then run DIME using something like:

```R
library(DIME)
hn1 = read.table("./data/vs_exprA_exprB_chr1")
hn1[is.na(hn1[,1]),1] = 0
pc = proc.time()
hn1.d1 = DIME(hn1[,1])
print(proc.time()-pc)
save(hn1.d1, file="./data/vs_exprA_exprB_500_chr1-diff")
#q(save="no")
```

But I must have messed up somewhere because even running a smaller chromosome I was hitting the max Wall time and not being able to finish running DIME

:(
