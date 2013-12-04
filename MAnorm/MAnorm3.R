# Author: Ying Wu daiyingw@gmail.com
# Last major update: December 2013
# This version is heavily modified from the original version
# to be faster and calculate p-value using edgeR
# License for this code should be the same as the original code but license is not given for original
# Latest version: https://github.com/ying-w/chipseq-compare/tree/master/MAnorm
# Original : http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/R_tutorial.html
#
# This script should be run with MAnorm3.sh in the same directory and is called 
# by this bash script at the 4th step
#
#####################################################################################################
# ideally we would just pass in a boolean vector for intersect instead of full matrix
# but because of the way that the files are generated in earlier step, this is easier
# prob will be more efficient with a complete rewrite
# Todo:
# finish fixing catch case for no replicates

# normalizeMA() code adapted from affy::normalize.loess function
# https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/affy/R/normalize.loess.R
normalizeMA = function(subset_mat, full_mat, method="loess", log.it=TRUE)
{
    # Normalize using MA line for each replicate
    # method can be either loess (NOT loWess) or rlm
    # should not return anything < 0 or else edgeR wont work

    #parameters from normalize.loess.R that I've modified/excluded
    # subsample=sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
    subsample = 1:nrow(subset_mat)
    epsilon = 10^-2
    maxit = 1
    verbose = TRUE

    J = ncol(full_mat) #dim(full_mat)[2]
    II = nrow(full_mat) #dim(full_mat)[1]
    if(log.it){
        subset_mat = log2(subset_mat)
        full_mat = log2(full_mat)
    }

    change = epsilon + 1
    iter = 0
    # weights for loess:
    # this way we give 0 weight to the
    # extremes added so that we can interpolate
    w = c(0, rep(1,length(subsample)), 0) 
    
    while(iter < maxit){
        iter = iter + 1
        means = matrix(0,II,J) ##contains temp of what we substract

        for (j in 1:(J-1)) {
            for (k in (j+1):J) {
                M = subset_mat[,j] - subset_mat[,k]
                A = (subset_mat[,j] + subset_mat[,k]) / 2
                ##put endpoints in so we can interpolate
                index = c(order(A)[1], subsample, order(-A)[1])
                M_subset = M[index]
                A_subset = A[index]
                if(method == "loess") {
                    span=2/3
                    family.loess="symmetric"
                    aux = loess(M_subset ~ A_subset, span=span, degree=1, weights=w, family=family.loess)
                } else if (method == "rlm") {
                    require(MASS)
                    aux = rlm(M_subset ~ A_subset)
                } else {
                    stop("Invalid method specified")
                }
                # plot(A_subset, M_subset)
                ## for rlm
                # abline(aux$coefficients)
                ## for loess
                # M_subset_rescale = predict(aux, data.frame(M_subset))
                # lines(A_subset, M_subset_rescale)
                ## plot changes
                # pp = predict(aux, data.frame(A_subset = A_subset))
                # plot(A_subset, M_subset); abline(aux$coefficients)#; lines(A_subset, pp)
                # plot(A_subset, M_subset-pp); abline(aux$coefficients); abline(h=0, col="blue")
                
                #colname of data frame below must match
                aux = predict(aux, data.frame(A_subset = (full_mat[,j] + full_mat[,k]) / 2)) / J
                means[, j] = means[, j] + aux
                means[, k] = means[, k] - aux
                if (verbose)
                message("Done with",j,"vs",k,"in iteration",iter)
            } #for k
        } #for j
        full_mat = full_mat - means
        change = max(colMeans((means[subsample,])^2))

        if(verbose)
        message("Iteration: ", iter, ", change: ", change)
    } #while

    if ((change > epsilon) & (maxit > 1)) warning(paste("No convergence after", maxit, "iterations."))

    if(log.it) {
        return(2^full_mat)
    } else
    return(full_mat)
}

extractTitle = function(titles, numsplit = 2) {
    # This function is to extract condition name from 
    # text field for title when plotting
    # possibilities are:
    #  title1, title2
    #  title1_rep1, title1_rep2, title2_rep1, title2_rep2
    #  title1, title2, title3 #with numsplit = 3
    # code written for readability not speed
    
    if(length(titles) %% numsplit) { stop(paste("Titles not a multiple of", numsplit)) }
    numreps = length(titles) / numsplit
    ret = rep("", numsplit)
    
    if(numreps == 1) {
        return(titles)
    }
    
    for(i in 1:numsplit) {
        current_title = titles[((i-1)*numreps+1):(i*numreps)]
        #make all of the titles the same length
        current_title_list = lapply(strsplit(current_title, ""), "[", 1:min(nchar(current_title)))
        current_subset = rep(FALSE, min(nchar(current_title)))
        #look for first mismatch since that is where you stop
        for(j in 2:numreps) { # numrep === length(current_title_list)
            current_subset = current_subset | current_title_list[[1]] != current_title_list[[j]]            
        }
        #current_subset will be TRUE if there is ANY mismatch
        #if NO mismatch then set end to end of match
        current_end = min(which(current_subset)-1, min(nchar(current_title)))
        if(current_end == 0) { ret[i] = paste0("Sample", i) }
        else { ret[i] = substr(current_title[1], 1,  current_end) }
    }
    ret
}

plotSmoothMA = function(mat, pvals = NULL, pval_cut = 0.01, plotfun = smoothScatter, 
    mergeReps = "append", sel=rep(FALSE, nrow(mat)), myTitle = NULL, prefix = "TF", 
    xlim = NULL, ylim = NULL, color = NULL, ...) {
    # Creates a smoothed MA plot
    # parts from pv.DBAplotMA() in DiffBind/R/analyze.R
    # mergeReps can be append or add, c("append", "add")
    # This will split the matrix into left and right with additional columns being replicates
    # use sel to obtain subset (say fdr cutoff)
    # takes UNlogged mat
    if(class(mat) == "DGELRT") { #edgeR glmLRT output
        # M = res$Fold
        # A = res$Conc
        # idx = res$FDR <= th
        M = mat$table$logFC
        A = mat$table$logCPM

        expr_name = extractTitle(rownames(mat$sample))
    }
    else if(class(mat) == "matrix" || class(mat) == "data.frame") {
        if(ncol(mat) > 2) {
            if(mergeReps == "append") { 
                message("Appending replicates")
                x = log2(mat[,c(1:(ncol(mat)/2))]+1)
                y = log2(mat[,c((ncol(mat)/2+1):ncol(mat))]+1)
                #only works since R is column order
                dim(x) = NULL
                dim(y) = NULL
            } else if(mergeReps == "add") {
                message("Adding replicates")
                x = log2(rowSums(mat[,c(1:(ncol(mat)/2))])+2)
                y = log2(rowSums(mat[,c((ncol(mat)/2):ncol(mat))])+2)
            } else { stop("Invalid mergeReps") }
        }
        M = x - y
        A = (x + y) / 2
        
        expr_name = extractTitle(colnames(mat))
        
    } else { stop("Unknown data type for mat, must be matrix or data frame") }

    # careful: sel calculation are scaled, if reps then sel is expanded to fill
    # length(A)/length(sel) is for scaling
    if(all(sel == FALSE) && !is.null(pvals)) {
        sel = pvals < pval_cut
        myTitle = sprintf('%s Binding Affinity: %s vs %s (%s of %s w/%s < %1.3f)',
            prefix, expr_name[1], expr_name[2], sum(sel)*(length(A)/length(sel)), length(A), "pval", pval_cut)
    } else if (sum(sel) > 0) {
        myTitle = sprintf('%s Binding Affinity: %s vs %s (%s of %s highlighted)',
            prefix, expr_name[1], expr_name[2], sum(sel)*(length(A)/length(sel)), length(A))
    } else {
        myTitle = sprintf('%s Binding Affinity: %s vs %s (%s sites)',
            prefix, expr_name[1], expr_name[2], length(A))
    }
    
    if(is.null(xlim)) { xlim = c(0, max(A)) }
    if(is.null(ylim)) { ylim = c(-max(abs(M)), max(abs(M))) }
    if(is.null(color) || length(color) == 1) { 
        cc = rep(rgb(0,0,0,0.3), length(A))
        if(is.null(color)) { color = rgb(1,0,0,0.5) }
        cc[sel] = color
    }
    
    #if length(A) is null then nothing will be displayed
    if(identical(plotfun, smoothScatter)) { #use identical() not ==
        plotfun(A, M, pch = 20, cex = 0.33, xlim = xlim, ylim = ylim,
            xlab = 'log concentration',
            ylab = sprintf('log2 fold change: %s - %s', expr_name[1], expr_name[2]),
            postPlotHook=grid(), #panel.first does not work here
            main = myTitle, panel.first=grid(), ...) 
        points(A[sel], M[sel], pch=20, cex=0.33, col=2) #highlight significant
    } else if(identical(plotfun, plot)) {
        plotfun(A, M, pch = 20, cex = 0.33, xlim = xlim, ylim = ylim,
            xlab = 'log concentration',
            ylab = sprintf('log2 fold change: %s - %s', expr_name[1], expr_name[2]),
            main = myTitle, panel.first=grid(), col = cc, ...)
    } else {
        message("Unknown plotting function, must use smoothScatter or plot\nPlotting skipped")
        return(invisible(cbind(M, A)))
    }

    abline(h=0,col='dodgerblue')
    
    if(sum(sel) > 0) {
        message(sum(sel)*(length(A)/length(sel)), "/", length(A), " or about ",
            round(sum(sel)*(length(A)/length(sel))/length(A)*100, 3), " % highlighted")
    }
    
    #In case you are wondering why everyone defines their own plotting function, see: http://xkcd.com/927/
    
    invisible(cbind(M, A))
}

#####################################################################################################
## BELOW IS CODE THAT IS RUN
#####################################################################################################
if(length(list.files(pattern="read1a.bed")) == 1) {  #11 arguments used, prob could be made more stringent
    # common is intersect, merge is pairwise overlaps merged
    common_peak_count_read1a = read.table("tmp_common_peak_read1a.counts",header=FALSE)
    common_peak_count_read2a = read.table("tmp_common_peak_read2a.counts",header=FALSE)
    common_peak_count_read1b = read.table("tmp_common_peak_read1b.counts",header=FALSE)
    common_peak_count_read2b = read.table("tmp_common_peak_read2b.counts",header=FALSE)
    peak_count_read1a = read.table("tmp_peak_read1a.counts",header=FALSE)
    peak_count_read2a = read.table("tmp_peak_read2a.counts",header=FALSE)
    peak_count_read1b = read.table("tmp_peak_read1b.counts",header=FALSE)
    peak_count_read2b = read.table("tmp_peak_read2b.counts",header=FALSE)
    merge_common_count_read1a = read.table("tmp_merge_common_read1a.counts", header=FALSE)
    merge_common_count_read1b = read.table("tmp_merge_common_read1b.counts", header=FALSE)
    merge_common_count_read2a = read.table("tmp_merge_common_read2a.counts", header=FALSE)
    merge_common_count_read2b = read.table("tmp_merge_common_read2b.counts", header=FALSE)
    merge_common_peak_count_read1a = read.table("tmp_merge_common_peak_read1a.counts",header=FALSE)
    merge_common_peak_count_read2a = read.table("tmp_merge_common_peak_read2a.counts",header=FALSE)
    merge_common_peak_count_read1b = read.table("tmp_merge_common_peak_read1b.counts",header=FALSE)
    merge_common_peak_count_read2b = read.table("tmp_merge_common_peak_read2b.counts",header=FALSE)
    common_count_mat = cbind(common_peak_count_read1a[,4], common_peak_count_read1b[,4], 
        common_peak_count_read2a[,4], common_peak_count_read2b[,4])
    all_count_mat = cbind(peak_count_read1a[,4], peak_count_read1b[,4], 
        peak_count_read2a[,4], peak_count_read2b[,4])
    common_merge_count_mat = cbind(merge_common_count_read1a[,4], merge_common_count_read1b[,4],
        merge_common_count_read2a[,4], merge_common_count_read2b[,4])
    all_merge_count_mat = cbind(merge_common_peak_count_read1a[,4], merge_common_peak_count_read1b[,4], 
        merge_common_peak_count_read2a[,4], merge_common_peak_count_read2b[,4])
} else { # no replicates
    common_peak_count_read1 = read.table("tmp_common_peak_read1.counts",header=FALSE)
    common_peak_count_read2 = read.table("tmp_common_peak_read2.counts",header=FALSE)
    peak_count_read1 = read.table("tmp_peak_read1.counts",header=FALSE)
    peak_count_read2 = read.table("tmp_peak_read2.counts",header=FALSE)
    merge_common_count_read1 = read.table("tmp_merge_common_read1.counts", header=FALSE)
    merge_common_count_read2 = read.table("tmp_merge_common_read2.counts", header=FALSE)
    merge_common_peak_count_read1 = read.table("tmp_merge_common_peak_read1.counts",header=FALSE)
    merge_common_peak_count_read2 = read.table("tmp_merge_common_peak_read2.counts",header=FALSE)
    common_count_mat = cbind(common_peak_count_read1, common_peak_count_read2)
    all_count_mat = cbind(peak_count_read1, peak_count_read2)
    common_merge_count_mat = cbind(merge_common=merge_common_count_read1, merge_common_count_read2)
    all_merge_count_mat = cbind(merge_common_peak_count_read1, merge_common_peak_count_read2)
}
table_MA = read.table("tmp_MAnorm.bed",header=FALSE)
table_merge_MA = read.table("tmp_MAnorm_merge.bed",header=FALSE)

#####################################################################################################
# new way + rescaling, works for 1 replicate only (0 replicates will give rowMeans error)
#####################################################################################################
#M = log2(rowMeans(common_count_mat[,1:(ncol(common_count_mat)/2)]+1)/
#    rowMeans(common_count_mat[,((ncol(common_count_mat)/2)+1):ncol(common_count_mat)]+1))
#M = log2((rowMeans(common_count_mat[,1:2])+1)/(rowMeans(common_count_mat[,3:4])+1)) #same as above
#A = rowMeans(log2(common_count_mat+1)) #same as below
#see rest of plotMA code
# pretty always returns a 2+5 vector
#smoothScatter(A,M,cex=1,main="MA plot before rescaling (common peaks)",  xlim=pretty(A)[c(1, 7)])

# name of files from name of folder (something like /protein_condition1_condition2)
f_path = strsplit(getwd(), "_")[[1]]
f_prefix = strsplit(getwd(), "/")[[1]]
f_prefix = f_prefix[length(f_prefix)]
if(length(f_path) > 2)
{
    x_name = paste0(rep(f_path[length(f_path)-1] ,ncol(all_count_mat)/2), c(1,2))
    y_name = paste0(rep(f_path[length(f_path)] ,ncol(all_count_mat)/2), c(1,2))
}
colnames(all_count_mat) = c(x_name, y_name)
colnames(common_count_mat) = c(x_name, y_name)

# MAnorm rescaling by subtracting MA adjustment
normalized_all_count_mat = normalizeMA(common_count_mat+1, all_count_mat+1, method="rlm")-1
#subset_mat = common_count_mat+1; full_mat = all_count_mat+1 #DEBUG

## plotting
#MA = plotSmoothMA(all_count_mat, plotfun="none")
#xlim = c(0, ceiling(max(MA[,2])))
#ylim = c(-ceiling(max(abs(MA[,1]))), ceiling(max(abs(MA[,1]))))

## MAplot before rescaling smoothed
#png(paste0(f_prefix, '_MAsmooth_before_rescaling.png'), width=3000, height=3000, res=300)
#plotSmoothMA(common_count_mat)
#dev.off()

## MAplot before rescaling as dots (not smoothed)
png(paste0(f_prefix, '_MAsmooth_before_rescaling_dots.png'), width=3000, height=3000, res=300)
plotSmoothMA(common_count_mat, plotfun=plot)
dev.off()

## If you want to see the MA line that is used:
png(paste0(f_prefix, '_MAsmooth_before_rescaling_dots_MAline.png'), width=3000, height=3000, res=300)
A_tmp = rowMeans(cbind(c(log2(all_count_mat[,1]+1), log2(all_count_mat[,2]+1)),
 c(log2(all_count_mat[,3]+1), log2(all_count_mat[,4]+1))))
offset_mat2 = log2(all_count_mat+1) - log2(normalized_all_count_mat+1)
M_tmp = c(offset_mat2[,1]-offset_mat2[,3], offset_mat2[,2]-offset_mat2[,4])
plotSmoothMA(common_count_mat, plotfun=plot)
points(A_tmp, M_tmp, pch=20, cex=0.2, col="blue")
dev.off()

## MAplot before rescaling (all peaks, above were common only)
#png(paste0(f_prefix, '_MAsmooth_before_rescaling_all.png'), width=3000, height=3000, res=300)
#plotSmoothMA(all_count_mat)
#dev.off()

# MAplot before rescaling w/common highlighted (all peaks)
png(paste0(f_prefix, '_MAsmooth_before_rescaling_all_common.png'), width=3000, height=3000, res=300)
plotSmoothMA(all_count_mat, sel = table_MA[,4] == "common_peak1" | table_MA[,4] == "common_peak2")
dev.off()

# MAplot before rescaling w/MA line (all peaks)
png(paste0(f_prefix, '_MAsmooth_before_rescaling_all_dots_MAline.png'), width=3000, height=3000, res=300)
plotSmoothMA(all_count_mat, plotfun=plot)
points(A_tmp, M_tmp, pch=20, cex=0.2, col="blue")
dev.off()

# Differential call using edgeR and normalized count matrix (subtract MA adjustment)
# DIFFBIND WAY (only works for replicates (group = 4))
# does not work as well as using offset matrix (see below)
if(0) { #old edgeR workflow
require(edgeR)
normalized_all_count_mat[normalized_all_count_mat < 0] = 0
res = DGEList(normalized_all_count_mat)
res$samples$group = c(0,0,1,1)
res$counts = round(res$counts) #method requires integers
res = calcNormFactors(res,method="TMM") #ignored when offset is specified
res$design = model.matrix(~res$samples$group)
res = estimateGLMCommonDisp(res,res$design)
# estimateGLMCommonDisp() calculates logCPM abundance and common dispersion
# logCPM is approx log2(rowMeans(cpm(resm))) (did not look up exact calculation)
# when this is run without estimateGLMTrendedDisp called, it will 
#  "squeeze tagwise dispersions towards common dispersion" edgeR 2.8.2
# this means that when using plotBCV it wont show the trend line
res = estimateGLMTagwiseDisp(res,res$design)  
res$GLM = glmFit(res,res$design)
res$LRT = glmLRT(res$GLM,2)
#out = topTags(res$LRT n = nrow(normalized_all_count_mat)) #res$LRT$table stores all the fun stuff
}

require(edgeR)
################################################################################
# Differential call using edgeR and normalized count matrix (offset matrix)
################################################################################
res = DGEList(all_count_mat, group = c(1,1,0,0))
offset_mat = log(all_count_mat+1) - log(normalized_all_count_mat+1)
# shift offset matrix: https://stat.ethz.ch/pipermail/bioconductor/2013-March/051680.html
# keep in mind that this ignores libsize diff since offset takes that into account
avgloglibsize = mean(log(res$samples$lib.size))
offset_mat = offset_mat - mean(offset_mat) + avgloglibsize 
#no TMM
res$design = model.matrix(~res$samples$group)
res = estimateGLMCommonDisp(res, res$design, offset = offset_mat)
res = estimateGLMTagwiseDisp(res, res$design, offset = offset_mat) #not sure if this should be included
#no trended
GLM = glmFit(res,res$design, offset = offset_mat); LRT = glmLRT(GLM,2)

# more plots
#plotSmoothMA(normalized_all_count_mat)
#png(paste0(f_prefix, '_MAsmooth_after_rescaling.png'), width=3000, height=3000, res=300)
#plotSmoothMA(LRT)
#dev.off()

require(qvalue) #possibly remove this and replace w/p.adjust
LRT$qv = qvalue(LRT$table$PValue)
summary(LRT$qv)
png(paste0(f_prefix, '_hist_pval.png'), width=3000, height=3000, res=300)
hist(LRT$table$PValue, breaks=100, main=paste("Pval distribution, pi0 =", 
    round(LRT$qv$pi0,3)), xlab="Pvalues calculated by edgeR")
dev.off()

#MAplots after rescale
png(paste0(f_prefix, '_MAsmooth_after_rescaling_sig.png'), width=3000, height=3000, res=300)
#plotSmoothMA(normalized_all_count_mat, pvals=LRT$qv$qvalues)
plotSmoothMA(LRT, pvals=LRT$qv$qvalues)
dev.off()

png(paste0(f_prefix, '_MAsmooth_after_rescaling_sig_dots.png'), width=3000, height=3000, res=300)
plotSmoothMA(LRT, pvals=LRT$qv$qvalues, plotfun=plot)
dev.off()

################################################################################
## do the same as above for merged dataset
################################################################################
colnames(common_merge_count_mat) = c(x_name, y_name)
colnames(all_merge_count_mat) = c(x_name, y_name)
normalized_all_merge_count_mat = normalizeMA(common_merge_count_mat+1, all_merge_count_mat+1, method="rlm")-1

#png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling.png'), width=3000, height=3000, res=300)
#plotSmoothMA(common_merge_count_mat)
#dev.off()

#png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling_dots.png'), width=3000, height=3000, res=300)
#plotSmoothMA(common_merge_count_mat, plotfun=plot)
#dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling_dots_MAline.png'), width=3000, height=3000, res=300)
A_tmp = rowMeans(cbind(c(log2(all_merge_count_mat[,1]+1), log2(all_merge_count_mat[,2]+1)),
 c(log2(all_merge_count_mat[,3]+1), log2(all_merge_count_mat[,4]+1))))
offset_mat2 = log2(all_merge_count_mat+1) - log2(normalized_all_merge_count_mat+1)
M_tmp = c(offset_mat2[,1]-offset_mat2[,3], offset_mat2[,2]-offset_mat2[,4])
plotSmoothMA(common_merge_count_mat, plotfun=plot)
points(A_tmp, M_tmp, pch=20, cex=0.2, col="blue")
dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling_all.png'), width=3000, height=3000, res=300)
plotSmoothMA(all_merge_count_mat)
dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling_all_common.png'), width=3000, height=3000, res=300)
plotSmoothMA(all_merge_count_mat, sel = table_merge_MA[,4] == "merged_common_peak")
dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling_all_dots_MAline.png'), width=3000, height=3000, res=300)
plotSmoothMA(all_merge_count_mat, plotfun=plot)
points(A_tmp, M_tmp, pch=20, cex=0.2, col="blue")
dev.off()

#Using offsets
resm = DGEList(all_merge_count_mat, group = c(1,1,0,0))
offset_mat2 = log(all_merge_count_mat+1) - log(normalized_all_merge_count_mat+1)
avgloglibsize = mean(log(resm$samples$lib.size))
offset_mat2 = offset_mat2 - mean(offset_mat2) + avgloglibsize
#no TMM
resm$design = model.matrix(~resm$samples$group)
resm = estimateGLMCommonDisp(resm,resm$design)
resm = estimateGLMTagwiseDisp(resm,resm$design)  
GLMm = glmFit(resm,resm$design, offset = offset_mat2)
LRTm = glmLRT(GLMm,2)

LRTm$qv = qvalue(LRTm$table$PValue)
summary(LRTm$qv)
png(paste0(f_prefix, '_hist_pval_merged.png'), width=3000, height=3000, res=300)
hist(LRTm$table$PValue, breaks=100, main=paste("Pval distribution (merged overlap), pi0 =", 
    round(LRTm$qv$pi0,3)), xlab="Pvalues calculated by edgeR")
dev.off()

#png(paste0(f_prefix, '_MAsmooth_merge_after_rescaling.png'), width=3000, height=3000, res=300)
#plotSmoothMA(LRTm)
#dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_after_rescaling_sig.png'), width=3000, height=3000, res=300)
plotSmoothMA(LRTm, pvals=LRTm$qv$qvalues)
dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_after_rescaling_sig_dots.png'), width=3000, height=3000, res=300)
plotSmoothMA(LRTm, pvals=LRTm$qv$qvalues, plotfun=plot)
dev.off()


if(0) { # old way for reference
#M = as.matrix(log2((common_peak_count_read1+1)/(common_peak_count_read2+1)))
#A = as.matrix(0.5*log2((common_peak_count_read1+1)*(common_peak_count_read2+1)))
# old way rescaling
#b=lm(M~A)$coefficients
#linear = lm(M~A)$coefficients
#b=robustRegBS(M,A,beta=linear)
b=rlm(M~A)$coefficients

png('MAplot_before_rescaling.png')
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
ma.plot(A,M,cex=1,main="MA plot before rescaling (common peaks)")
abline(b,col="green")
dev.off()

cat("M = b[1] + b[2] * A\n")
log2_peak_count_read1 = log2(peak_count_read1 + 1)
log2_peak_count_read2 = log2(peak_count_read2 + 1)
log2_peak_count_read1_rescaled = (2-b[2])*log2_peak_count_read1/(2+b[2]) - 2*b[1]/(2+b[2]);
#you get the above forumla after solving for read1 when you plug in M and A
M_rescaled = (log2_peak_count_read1_rescaled - log2_peak_count_read2);
A_rescaled = (log2_peak_count_read1_rescaled + log2_peak_count_read2)/2;

png('MAplot_after_rescaling.png')
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
ma.plot(as.matrix(A_rescaled),as.matrix(M_rescaled),cex=1,main=" MA plot after rescaling (all peaks)")
dev.off ()

log2_merge_common_peak_count_read1 = log2(merge_common_peak_count_read1 + 1)
log2_merge_common_peak_count_read2 = log2(merge_common_peak_count_read2 + 1)
log2_merge_common_peak_count_read1_rescaled = (2-b[2])*log2_merge_common_peak_count_read1/(2+b[2]) - 2*b[1]/(2+b[2]);
merge_M_rescaled = (log2_merge_common_peak_count_read1_rescaled - log2_merge_common_peak_count_read2);
merge_A_rescaled = (log2_merge_common_peak_count_read1_rescaled + log2_merge_common_peak_count_read2)/2;
}

#output (only works for replicates)

table_MA[,5] = peak_count_read1a[,4]
table_MA[,6] = peak_count_read1b[,4]
table_MA[,7] = peak_count_read2a[,4]
table_MA[,8] = peak_count_read2b[,4]
table_MA[,9] = res$counts[,1]
table_MA[,10] = res$counts[,2]
table_MA[,11] = res$counts[,3]
table_MA[,12] = res$counts[,4]
table_MA[,13] = LRT$table$logFC
table_MA[,14] = LRT$table$logCPM
table_MA[,15] = -log10(LRT$table$PValue)
table_MA[,16] = -log10(LRT$qv$qvalues)

colnames(table_MA)[1] = "chr"
colnames(table_MA)[2] = "start"
colnames(table_MA)[3] = "end"
colnames(table_MA)[4] = "description"
colnames(table_MA)[5] = "#raw_read_1a"
colnames(table_MA)[6] = "#raw_read_1b"
colnames(table_MA)[7] = "#raw_read_2a"
colnames(table_MA)[8] = "#raw_read_2b"
colnames(table_MA)[9] = "#MAnorm_read_1a"
colnames(table_MA)[10] = "#MAnorm_read_1b"
colnames(table_MA)[11] = "#MAnorm_read_2a"
colnames(table_MA)[12] = "#MAnorm_read_2b"
colnames(table_MA)[13] = "logFC"
colnames(table_MA)[14] = "logCPM"
colnames(table_MA)[15] = "-log10(p-value)"
colnames(table_MA)[16] = "-log10(q-value)"

write.table(table_MA, paste0(f_prefix,"_MAnorm_result.xls"), sep="\t", quote=FALSE, row.names=FALSE)
#can adjusting make numnbers less than 0?
#adjusting can boost 0s to non-0
#careful about using qvalues package because diff pi0 can affect adjustment
#look at distribution of pvalues

# table_merge
table_merge_MA[,5] = merge_common_peak_count_read1a[,4]
table_merge_MA[,6] = merge_common_peak_count_read1b[,4]
table_merge_MA[,7] = merge_common_peak_count_read2a[,4]
table_merge_MA[,8] = merge_common_peak_count_read2b[,4]
table_merge_MA[,9] = resm$counts[,1]
table_merge_MA[,10] = resm$counts[,2]
table_merge_MA[,11] = resm$counts[,3]
table_merge_MA[,12] = resm$counts[,4]
table_merge_MA[,13] = LRTm$table$logFC
table_merge_MA[,14] = LRTm$table$logCPM
table_merge_MA[,15] = -log10(LRTm$table$PValue)
table_merge_MA[,16] = -log10(LRTm$qv$qvalues)

colnames(table_merge_MA)[1] = "chr"
colnames(table_merge_MA)[2] = "start"
colnames(table_merge_MA)[3] = "end"
colnames(table_merge_MA)[4] = "description"
colnames(table_merge_MA)[5] = "#raw_read_1a"
colnames(table_merge_MA)[6] = "#raw_read_1b"
colnames(table_merge_MA)[7] = "#raw_read_2a"
colnames(table_merge_MA)[8] = "#raw_read_2b"
colnames(table_merge_MA)[9] = "#MAnorm_read_1a"
colnames(table_merge_MA)[10] = "#MAnorm_read_1b"
colnames(table_merge_MA)[11] = "#MAnorm_read_2a"
colnames(table_merge_MA)[12] = "#MAnorm_read_2b"
colnames(table_merge_MA)[13] = "logFC"
colnames(table_merge_MA)[14] = "logCPM"
colnames(table_merge_MA)[15] = "-log10(p-value)"
colnames(table_merge_MA)[16] = "-log10(q-value)"

write.table(table_merge_MA, paste0(f_prefix,"_MAnorm_merged_result.xls"), sep="\t", quote=FALSE, row.names=FALSE)

rm(table_MA, table_merge_MA) #Using R CMD BATCH saves workspace so remove redundant tables
# source("../edgeRmatricies.R") #for additional figures for comparison of normalization
