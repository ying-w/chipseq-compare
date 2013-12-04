# DIME

Here is some code to preprocess data for [DIME](http://www.stat.osu.edu/~statgen/SOFTWARE/DIME/).
Their 2011 Bioinformatics paper can be found [here](http://pubmed.gov/21471015).

The author (C. Taslim) was kind enough to send me some matlab code to do the loess normalization used as input for DIME.
The matlab code can be found in the following 3 files and implements normalization from their 2009 Bioinformatics paper [here](http://pubmed.gov/19561022)

    normalizeMeanVarStep1.m
    normalizeMeanVarStep2.m
    runNormzMeanVarStep1_2.m
    
I have rewritten parts the three matlab scripts into an R file: `loess.R`
These 3 matlab files will not run under [Octave](http://www.gnu.org/software/octave/) due the `smooth()` function missing a loess option.
The output of this script is used as input for DIME.

If memory serves me correctly, I then made non-overlapping bins over an entire chromosome and counted reads in each bin using something like:

```bash
perl covBEDforDIME.pl chromsizes19M.tab > 500bp_window.bed &

intersectBed -v -abam ../experimentA1.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprA1.cov 
intersectBed -v -abam ../experimentA2.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprA2.cov 
intersectBed -v -abam ../experimentB1.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprB1.cov
intersectBed -v -abam ../experimentB2.bam -b ../artifact2.bed | coverageBed -abam stdin -b 500bp_window.bed | sortBed > ./data/exprB2.cov
```
`chromsizes19M.tab` is a two colum file of chromosome names in the first column and chromosome sizes in the second column
 
`artifact2.bed` is merged file containing DAC Blacklisted Regions and Duke Excluded Regions from ENCODE project

I then split up the file the coverage files by chromosome using `splitcov.pl` and run `loess.R` and then run DIME (below):

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

Here is a couple of lines from the `vs_exprA_exprB_chr1` file, this file is 118748 lines long:

    1.33995979450539
    1.23717107628172
    0
    -1.33995979450539
    3.50775704404255
    1.90947995192546
    0
    0

## Future direction
Probably not going to work on this further but one could cut out certain regions (peaks) 
or make the window size larger or take half a chromosome but for now I will investigate other techniques.
