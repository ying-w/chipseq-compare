# Comparing edgeR, DiffBind and DBChIP

Comparisons between there three programs are not straightforward and I will 
highlight here some key differences and code to get around these issues.
I have a paper coming out with more details about how these methods differ.

# Input format
Typically at this step, two types of files are required:

1.	Reads (bam/bed)
2.	Peaks (bed)

Peaks define the regions of interest which are then counted over in reads

## edgeR - [website](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
This program has the simplest input. It requires a counts matrix with features 
(peaks/genes) as rows and samples as columns. Replicates are necessary to 
estimate dispersion, however, dispersion could also be user supplied. 
Additionally, a design matrix must be specified by the user for this method.

The counts matrix can be created using counting over bam files using 
`intersectBed` from [bedtools]( http://bedtools.readthedocs.org/en/latest/), 
`htseq-count` from [htseq](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html),
`summarizeOverlaps` from [GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html) 
(and soon will be moved into [GenomicAlignments]( http://www.bioconductor.org/packages/devel/bioc/html/GenomicAlignments.html)),
and `featureCount` from [Rsubread]( http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html).
All these methods are similar in concept but may have slight differences in practice due to different implementations handling edge cases differently.
The last two are R packages and `featureCount` is supposed to be more [memory efficient](http://www.ncbi.nlm.nih.gov/pubmed/24227677).

## DiffBind - [website](http://www.bioconductor.org/packages/release/bioc/html/DiffBind.html)
This program requires csv file (sample sheet) that specifies your samples and their locations. 
The format of this csv file is detailed in the manual but a quick overview is as follows:

`SampleID,Tissue,Factor,Condition,Replicate,bamReads,bamControl,Peaks,PeakCaller`

The first four columns are attributes for your samples that will be used to create a design matrix while the last four columns are for specifying file locations and type. The read files can be gzipped bed files or bam files.

Design matrix will be generated (as contrast) based on your attributes, based on the way that DiffBind is currently coded it is not possible to insert your own design matrix. It is possible to do multi-factor comparisons using `block` in `dba.contrast()` but seems limited in functionality. When doing pairwise analysis, fold changes might be the opposite of what you expect, this occurs due to alphanumeric sorting at one step (somewhere in the code `levels()` is done and that’s the order for comparison) so to avoid this, add number in front of the name of the condition when specified during dba.contrast(), higher number will be subtracted from the lower number when sorted alphanumerically.

Counting will be done in parallel by custom C++ code (croi) or `summarizeOverlaps` (for low memory).
The C++ code will also extend reads based on `insertLength` and compute library size.

```r
library(DiffBind)
DBA = dba(sampleSheet="myexperiment.csv") 
Sys.sleep(1) #sometimes you will need this when parallel reading
DBA = dba.count(DBA, bCorPlot=FALSE) #surpress plots when on server
DBA = dba.contrast(DBA, group1=c(1,2), group2=c(3,4), name1="1cond", name2="2cond") 
# group1 and group2 referring to rows from myexperiment.csv
# in this case, I have two conditions with two replicates each
# fold changes will be for 2cond-1cond 
```

## DBChIP - [website](http://www.bioconductor.org/packages/release/bioc/html/DBChIP.html)
The code for this program is quite short and easy to follow, there is a DBChIP.R which is all the functions except NCIS and formats the data in such a way that the functions in NCIS.R will work.

There are two practical ways to import your data into DBChIP. Either format it as something they call MCS (minimum chip seq, which is pretty much a BED file) or use the [ShortRead]( http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html) to read in your data and feed in an AlignedRead object. The last format is to specify a BED file but I found that inconvenient because I do not keep around non-compressed BED files. I opted for using MCS format (below).

```r
library(DBChIP)
rawreads = list()

tmp = read.table("cond1rep1.bed.gz")
fiveend = tmp[,2] #take start
fiveend[tmp[,6] == "-"] = tmp[tmp[,6] == "-",3] #replace w/end for minus strand
# wish I could figure out a way to do this in a one line apply statement
rawreads[["cond1_r1"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)

tmp = read.table("cond1rep2.bed.gz")
fiveend = tmp[,2] #take start
fiveend[tmp[,6] == "-"] = tmp[tmp[,6] == "-",3] 
rawreads[["cond1_r2"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)

tmp = read.table("cond2rep1.bed.gz")
fiveend = tmp[,2] #take start
fiveend[tmp[,6] == "-"] = tmp[tmp[,6] == "-",3] 
rawreads[["cond2_r1"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)

tmp = read.table("cond2rep2.bed.gz")
fiveend = tmp[,2] #take start
fiveend[tmp[,6] == "-"] = tmp[tmp[,6] == "-",3] 
rawreads[["cond2_r2"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)

inputreads = list() #reads from background/reverse crosslink/input 

#names(inputreads) must match names(rawreads) for older versions of DBChIP
#with newer versions of DBChIP you can avoid this redundancy by specifying matching.input.names
tmp = read.table("cond1input.bed.gz")
fiveend = tmp[,2] #take start
fiveend[tmp[,6] == "-"] = tmp[tmp[,6] == "-",3] 
inputreads[["cond1_r1"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)
inputreads[["cond1_r2"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)

tmp = read.table("cond2input.bed.gz")
fiveend = tmp[,2] #take start
fiveend[tmp[,6] == "-"] = tmp[tmp[,6] == "-",3] 
inputreads[["cond2_r1"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)
inputreads[["cond2_r2"]] = data.frame(chr=paste(tmp[,1]), strand= tmp[,6], pos=fiveend)

# load peaks
# col7 is the signal column for peak caller output (bigger number = bigger peak) 
# col10 is the summit offset from peak caller output
# pos will hold position of summit (start + offset)
binding.site.list = list()
tmp = read.table("cond1peaks.bed”)
binding.site.list[["cond1"]] = data.frame(chr=paste(tmp[,1]), pos=tmp[,2]+tmp[,10], weight=tmp[,7])
tmp = read.table("cond2peaks.bed”)
binding.site.list[["cond2"]] = data.frame(chr=paste(tmp[,1]), pos=tmp[,2]+tmp[,10], weight=tmp[,7])
bs.list = read.binding.site.list(binding.site.list)
consensus.site = site.merge(bs.list) #merging sites from different conditions

conds = factor(c("cond1","cond1","cond2","cond2")) #cannot use more than 2 at the same time, gives error at test.diff.binding() step
cond1_vs_cond2_dat = load.data(chip.data.list=chip.data.list, conds=conds, 
    consensus.site=NULL, input.data.list=inputreads, #matching.input.names=matching.input.names,
    data.type="MCS", frag.len=round(mean(c(readlength[["cond1"]], readlength[["cond2"]]))))
# readlength is a list that holds fragment length
# the frag.len variable is used as shift.len (half of frag.len) to shift the windows when scanning the genome using windows
cond1_vs_cond2_dbchip$consensus.site = consensus.site #force NCIS to be used in load.data
cond1_vs_cond2_dbchip = get.site.count(cond1_vs_cond2_dbchip)
```

In this case, I used the same files as I specified in myexperiment.csv for DiffBind. Even though DiffBind will import the reads and shifts them, it does not save the read data so you cannot just load it into DBChIP. DBChIP tends to reduce the amount of information used as input (only taking 5’ end in MCS format for reads, only taking summit location for peaks).


# Reference binding regions used
This section will go over how peaks that are provided are merged / combined to define reference regions that are counted over. Ie. How to select features to count over since peaks in different conditions will overlap.

## edgeR
Since counts matrix is user specified, the user would also have to make the decision about what they want to use for reference binding regions. 

## DiffBind
All peaks provided in the csv file will be merged. There is a `minOverlap` tuning parameter. There is no option for user defined reference to count over. The merging will be done during sample sheet readin (`dba()`) using the private function `pv.peakset()`.

## DBChIP
User provides peak summits (one coordinate not start/end) and these summits will be combined clustering (`site.merge()`). Default parameters merge centroids within 100bp and keeps separate centroids greater than 250bp away. Do not specify `consensus.site` in `load.data()` if you want to use NCIS. Note that output for many DBChIP functions are lists by chromosome which can be annoying sometimes (but easy to `lapply()` over).

# Evaluation options
## edgeR
How to run edgeR is pretty well documented in the [user guide](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf). Briefly, the usual workflows are:

    # pairwise
    y <- DGEList(counts=x,group=factor(c(1,1,2,2)))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    topTags(et)

or

    #GLM
    design <- model.matrix(~factor(c(1,1,2,2)))
    y <- estimateGLMCommonDisp(y,design)
    y <- estimateGLMTrendedDisp(y,design)
    y <- estimateGLMTagwiseDisp(y,design)
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    topTags(lrt)

Trended Dispersion estimate might not be applicable for ChIP-seq and it would be good idea to check fits using `plotBCV()`

## DiffBind
This program has 3 different evaluation methods (edgeR/DESeq/DESeq2). By default, edgeR is run in parallel on the contrasts specified (code in previous step) using tagwise dispersion. (See Technical Notes part of DiffBind User Guide for more details)

```r
cond1_vs_cond2_diffbind  = dba.analyze(DBA, bCorPlot=FALSE)
head(dba.report(cond1_vs_cond2_diffbind )) #show top 6 as genomic ranges
```

## DBChIP
Similar to DiffBind (above), the command to do differential testing is only one line. A variance estimation method will be performed for conditions with no replicates (see their paper for details).  

```r
cond1_vs_cond2_dbchip = test.diff.binding(cond1_vs_cond2_dbchip)
```

# Comparing
Since all 3 of the programs covered use edgeR, you would think comparing outputs are simple; however, differences in reference binding regions will manifest as different number of rows for binding matrix. To get around this, you could modify the starting peaks.
 
## edgeR
One could use the binding matrix from DiffBind or DBChIP though it is kind of pointless in the case of DiffBind because the `bReduceObjects` option in `dba.analyze()` would return the edgeR object. The edgeR analysis is what you would get without subtracting input. Default uses effective library size normalization.

## DiffBind
Input is scaled (down only) and subtracted from counts. Default is to use full library size normalization.

    scale_factor = colSums(counts) / colSums(counts_input) #scale = cond$libsize / cont$libsize
    for(i in length(scale_factor)) {
        if(scale_factor[i] > 1) scale_factor[i] = 1 # dont upscale control
        if(scale_factor[i] != 0) {
           counts_input[,i] = ceiling(counts_input[,i] * scale_factor[i])
        }
    }

## DBChIP
Input can be scaled using NCIS before subtraction. A median ratio is calculated and used as library size for normalization.

    # taken from DBChIP::median.ratio()
    gm <- exp(rowMeans(log(counts)))
    return(apply(counts, 2, function(x) median((x/gm)[gm > 0])))


    # something like this
    scalefactorNCIS = function(tmplist) {
        sapply(1:ncol(tmplist$site.count), function(x) DBChIP:::NCIS.internal(
            tmplist$chip.list[[x]], tmplist$input.list[[x]], 
            chr.vec=tmplist$chr.vec, chr.len.vec=tmplist$chr.len.vec)$est)
    }

# Code for subtracting and normalizing differently using NCIS output
Choices are as follows:
- library size (full/median/effective)
- Input scaling and subtraction (DiffBind/NCIS/none)
- Model input as offsets instead of subtracting from counts (yes/no)

DiffBind default is : full/DiffBind/no (for newer versions) and effective/DiffBind/no (for older versions)
DBChIP default is : median/NCIS/no when NCIS is run

`libsize` will need to be defined

    get.library.size = function(counts, f_name, method="effective") {
        # method can be full/median/effective
        message(paste("libnorm =", method))
        
        libsize = list()
        libsize[["gr_high"]] = c(19291260, 16754796)
        
        ## full library size (DiffBind)
        if(method == "full") { return(c(libsize[[f_name[1]]], libsize[[f_name[2]]])) }
        
        ## median library size (DBChIP)
        if(method == "median") { # taken from DBChIP::median.ratio()
            gm <- exp(rowMeans(log(counts)))
            return(apply(counts, 2, function(x) median((x/gm)[gm > 0])))
        }

        ## effective library size
        if(method == "effective") { return(colSums(counts)) }
        
    }

    scalefactorNCIS = function(tmplist) {
        sapply(1:ncol(tmplist$site.count), function(x) DBChIP:::NCIS.internal(
            tmplist$chip.list[[x]], tmplist$input.list[[x]], 
            chr.vec=tmplist$chr.vec, chr.len.vec=tmplist$chr.len.vec)$est)
    }
    # input scaling and subtraction
    get.scaling.factor = function(counts, counts_input, scaleinput="nosub", return.final = TRUE, dat=NULL, chipBED=NULL, NCISscale=NULL) {
    # scaleinput can be DiffBind/NCIS/nosub
    # digging around https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/DiffBind/ (rev 80867)

        ## subtract input after scaling (DiffBind)
        if(scaleinput == "DiffBind") { 
        # from counts.R:pv.counts()
        # libsize is created by DiffBind:::pv.getCounts which calls out to C for call to Ptr->size()
        # guess that libsize is just sum
        # cont = control, cond = condition
        # edgeR call is in analyze.R:pv.DEinit() and calls utils.R:pv.get_reads()
         
            message("DiffBind scaling")
            # from pv.counts()
            scale_factor = colSums(counts) / colSums(counts_input) #scale = cond$libsize / cont$libsize
            for(i in length(scale_factor)) {
                if(scale_factor[i] > 1) scale_factor[i] = 1 # dont upscale control
                if(scale_factor[i] != 0) {
                   counts_input[,i] = ceiling(counts_input[,i] * scale_factor[i])
                }
            }
            # summarized from pv.get_reads()
            counts_final = pmax(counts - counts_input, 1)
        } #DiffBind
        
        ## subtract input after NCIS (DBChIP)
        if(scaleinput == "NCIS") { 
            # from the following functions in DBChIP
            # get.site.count() and get.NCIS.norm.factor()
            
            message("NCIS scaling")
            # BED reading is not tested
            # in code, it is read in as a list since multiple conditions then using lapply on read.MCS
            # must be in chr/strand/pos
            if(!is.null(dat)) {
                tmplist = dat
            } else if(is.null(dat) && !is.null(chipBED)) {
                fiveend = chipBED[,2]; fiveend[chipBED[,6] == "-"] = chipBED[chipBED[,6] == "-",3]
                tmpchip = DBChIP:::read.MCS(data.frame(chr=paste(chipBED[,1]), strand=chipBED[,6], pos=fiveend))
                fiveend = inputBED[,2]; fiveend[inputBED[,6] == "-"] = inputBED[inputBED[,6] == "-",3]
                tmplist = DBChIP:::read.MCS(data.frame(chr=paste(inputBED[,1]), strand=inputBED[,6], pos=fiveend))
            }
            
            if(!is.null(NCISscale)) {
                scale_factor = NCISscale
            } else { 
                scale_factor = scalefactorNCIS(tmplist) 
            } # own function since calculation take a while so allow for precomputation
            
            counts_input = t(t(counts_input)*scale_factor)
            counts_final = counts - counts_input
            counts_final = pmax(round(counts_final), 0)
        }

        ## no subtraction
        if(scaleinput == "nosub" || scaleinput == "edgeR") { message("No input scaling"); counts_input = NULL; counts_final = counts }

        if(return.final) { counts_final } else { counts_input }
    }

    calcLRT = function(counts_final, cond1, cond2, libnorm, counts_input = NULL) {
    # libnorm can be full/median/effective

        if(is.null(counts_input)) {
            message("No offset used")
            z = DGEList(counts_final+1, group = c(1,1,0,0), lib.size = get.library.size(counts_final+1, c(cond1, cond2), method=libnorm))
            #no TMM
            z$design = model.matrix(~z$samples$group)
            z = estimateGLMCommonDisp(z, z$design)
            z = estimateGLMTagwiseDisp(z, z$design)
            #no trended
            LRT = glmLRT(glmFit(z,z$design),2)
        } else {
            message("Offset used") # input as offset
            z = DGEList(counts_final+1, group = c(1,1,0,0), lib.size = get.library.size(counts_final+1, c(cond1, cond2), method=libnorm))
            offset_mat = t(matrix(rep(log(z$samples$lib.size),nrow(counts_input)),nrow=ncol(counts_input))) + log(counts_input+1)
            #no TMM
            z$design = model.matrix(~z$samples$group)
            z = estimateGLMCommonDisp(z, z$design, offset = offset_mat)
            z = estimateGLMTagwiseDisp(z, z$design, offset = offset_mat) #not sure if offset should be reincluded
            #no trended
            LRT = glmLRT(glmFit(z,z$design, offset = offset_mat),2)
        }
    }

    getInputCounts = function(dat, window.size=250) {
    # modified from DBChIP::get.site.count
        consensus.site <- dat$consensus.site
        res = matrix(0, nrow = sum(unlist(lapply(dat$consensus.site, nrow))), ncol = length(names(dat$chip.list)))
        currow = 1
        for(chr in names(consensus.site)){
            hi <- DBChIP:::get.site.count.hist(consensus.site[[chr]]$pos, window.size=window.size)
            #chip.count.vec <- sapply(dat$chip.list, function(x) get.chr.site.count(x[[chr]], hi))
            
            if(length(consensus.site[[chr]]$pos)==1) chip.count.vec <- t(as.matrix(chip.count.vec))
            #if(!is.null(dat$input.list) && subtract.input){
                input.count.vec <- sapply(names(dat$chip.list), function(x) DBChIP:::get.chr.site.count(dat$input.list[[dat$matching.input.names[x]]][[chr]], hi))
            #    chip.count.vec <- chip.count.vec - t(t(input.count.vec)*dat$norm.factor.vec)
            #    chip.count.vec <- pmax(round(chip.count.vec), 0)
            #}
            #rownames(input.count.vec) <- paste(chr, "_", consensus.site[[chr]]$pos, sep="")
            #res <- rbind(res, chip.count.vec)
            res[currow:(currow+nrow(input.count.vec)-1),] = input.count.vec #would be more efficient to do as col and t()
            currow = currow+nrow(input.count.vec)
        }
        #dat$site.count <- res
        return(res)
    }

    getDiff = function(dat, libnorm, f_prefix, scaleinput, useOffset=FALSE, returnLRT = FALSE, NCISscale=NULL) {
        # libnorm can be full/median/effective
        # scaleinput can be DiffBind/NCIS/nosub
        
        f_split = strsplit(f_prefix, "_")[[1]]
        counts = dat$site.count
        counts_input = getInputCounts(dat) # subtract raw counts
        if(useOffset) {
            count_adjust = get.scaling.factor(counts, counts_input, scaleinput=scaleinput, return.final = FALSE, dat = dat, NCISscale=NCISscale)
            LRT = calcLRT(counts, paste(f_split[1], f_split[2], sep = "_"), paste(f_split[1], f_split[3], sep = "_"), libnorm, count_adjust)
        } else { #offsets not used
            counts_final = get.scaling.factor(counts, counts_input, scaleinput=scaleinput, return.final = TRUE, dat = dat, NCISscale=NCISscale)
            LRT = calcLRT(counts_final, paste(f_split[1], f_split[2], sep = "_"), paste(f_split[1], f_split[3], sep = "_"), libnorm)
        }
        DTD = decideTestsDGE(LRT, p.value=pvalcut)
        message(sum(DTD != 0))
        if(returnLRT) {
            invisible(LRT)
        } else {
            invisible(sum(DTD != 0))
        }
    }

    runallmethods = function(dat, f_prefix) {
        res = data.frame("DiffBind" = rep(0,6), "NCIS" = rep(0,6), "nosub" = rep(0,6), 
            row.names=c("full-offset", "full-nooffset", "median-offset", "median-nooffset", "effective-offset", "effective-nooffset"))
        NCISscale = scalefactorNCIS(dat)
        for(scaleinput in c("DiffBind", "NCIS", "nosub")) 
            for(libnorm in c("full-offset", "full-nooffset", "median-offset", "median-nooffset", "effective-offset", "effective-nooffset"))
                if(strsplit(libnorm, "-")[[1]][2] == "offset") {
                    if(scaleinput != "nosub")
                        res[libnorm, scaleinput] = getDiff(dat = dat, libnorm = strsplit(libnorm, "-")[[1]][1], f_prefix = f_prefix, scaleinput = scaleinput, useOffset = TRUE, NCISscale=NCISscale) 
                } else {
                    res[libnorm, scaleinput] = getDiff(dat = dat, libnorm = strsplit(libnorm, "-")[[1]][1], f_prefix = f_prefix, scaleinput = scaleinput, useOffset = FALSE, NCISscale=NCISscale) 
                }
        res
    }



    > runallmethods(gr_hl_dbchip, f_prefix = "gr_high_low")
    #                    DiffBind  NCIS nosub
    # full-offset           17222 17237     0
    # full-nooffset         17240 17239 17222
    # median-offset          3287  2614     0
    # median-nooffset        2943  3597  3287
    # effective-offset       3538  2776     0
    # effective-nooffset     2983  3392  3538

# Appendix
```
> sessionInfo()

## R version 2.15.3 (2013-03-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=C                 LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] Vennerable_2.2       xtable_1.7-1         gtools_3.0.0        
##  [4] reshape_0.8.4        plyr_1.8             RColorBrewer_1.0-5  
##  [7] RBGL_1.34.0          graph_1.36.2         DBChIP_1.2.0        
## [10] DESeq_1.10.1         lattice_0.20-23      locfit_1.5-9.1      
## [13] DiffBind_1.4.2       Biobase_2.18.0       GenomicRanges_1.10.7
## [16] IRanges_1.16.6       BiocGenerics_0.4.0   edgeR_3.0.8         
## [19] limma_3.14.4         knitr_1.4.1         
## 
## loaded via a namespace (and not attached):
##  [1] amap_0.8-7           annotate_1.36.0      AnnotationDbi_1.20.7
##  [4] codetools_0.2-8      DBI_0.2-7            digest_0.6.3        
##  [7] evaluate_0.4.7       formatR_0.9          gdata_2.13.2        
## [10] genefilter_1.40.0    geneplotter_1.36.0   gplots_2.11.0       
## [13] RSQLite_0.11.4       splines_2.15.3       stats4_2.15.3       
## [16] stringr_0.6.2        survival_2.37-4      tools_2.15.3        
## [19] XML_3.98-1.1         zlibbioc_1.4.0
```
