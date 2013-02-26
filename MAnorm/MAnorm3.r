# major modifications to MAnorm to be faster + calculate p-value using edgeR
#####################################################################################################
# ideally we would just pass in a boolean subset index instead of full matrix
# but because of the way that the files are generated in earlier step, this is easier
# code mostly copied from affy::normalize.loess function
# https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/affy/R/normalize.loess.R
normalizeMA = function(subset_mat, full_mat, method="loess", log.it=TRUE)
{
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
                if(method == "loewss") {
                    span=2/3
                    family.loess="symmetric"
                    aux = loess(M_subset ~ A_subset, span=span, degree=1, weights=w, family=family.loess)
                } else if (method == "rlm") {
                    require(MASS)
                    aux = rlm(M_subset ~ A_subset)
                } else {
                    stop("Invalid method specified")
                }
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

# parts from pv.DBAplotMA() in DiffBind/R/analyze.R
# mergeReps can be append or add, c("append", "add")
# use sel to obtain subset (say fdr cutoff)
plotSmoothMA = function(mat, pvals = NULL, pval_cut = 0.01, plotfun = smoothScatter, 
    mergeReps = "append", sel=rep(FALSE, nrow(mat)), myTitle = NULL) {
    #takes UNlogged mat
    if(class(mat) == "DGELRT") {
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
        # M = res$Fold
        # A = res$Conc
        # idx = res$FDR <= th
        M = x - y
        A = (x + y) / 2
        
        expr_name = extractTitle(colnames(mat))
        
    } else { stop("Unknown matrix type") }

    # todo: sel calculation below is wrong, if reps then sel is expanded to fill
    if(all(sel == FALSE) && !is.null(pvals)) {
        sel = pvals < pval_cut
        myTitle = sprintf('%s Binding Affinity: %s vs %s (%s of %s w/%s < %1.3f)',
            'TF', expr_name[1], expr_name[2], sum(sel), length(A), "pval", pval_cut)
    } else if (sum(sel) > 0) {
        myTitle = sprintf('%s Binding Affinity: %s vs %s (%s of %s highlighted)',
            'TF', expr_name[1], expr_name[2], sum(sel), length(A))
    } else {
        myTitle = sprintf('%s Binding Affinity: %s vs %s (%s sites)',
            'TF', expr_name[1], expr_name[2], length(A))
    }
    if(identical(plotfun, smoothScatter)) { #use identical() not ==
        plotfun(A, M, pch = 20, cex = 0.33, xlim = c(0, ceiling(max(A))),
            xlab = 'log concentration',
            ylab = sprintf('log fold change: %s - %s', expr_name[1], expr_name[2]),
            postPlotHook=grid(), 
            main = myTitle, panel.first=grid()) #if length(A) is null then nothing will be printed
    } else {
        plotfun(A, M, pch = 20, cex = 0.33, xlim = c(0, ceiling(max(A))),
            xlab = 'log concentration',
            ylab = sprintf('log fold change: %s - %s', expr_name[1], expr_name[2]),
            main = myTitle, panel.first=grid()) #if length(A) is null then nothing will be printed
    }

    abline(h=0,col='dodgerblue')
    points(A[sel], M[sel], pch=20, cex=0.33, col=2) #highlight significant
}

#####################################################################################################
## Run this code
#####################################################################################################
if(length(list.files(pattern="read1a.bed")) == 1) {  #11 arguments used, prob could be made more stringent
    common_peak_count_read1a = read.table("tmp_common_peak_count_read1a",header=FALSE)
    common_peak_count_read2a = read.table("tmp_common_peak_count_read2a",header=FALSE)
    common_peak_count_read1b = read.table("tmp_common_peak_count_read1b",header=FALSE)
    common_peak_count_read2b = read.table("tmp_common_peak_count_read2b",header=FALSE)
    peak_count_read1a = read.table("tmp_peak_count_read1a",header=FALSE)
    peak_count_read2a = read.table("tmp_peak_count_read2a",header=FALSE)
    peak_count_read1b = read.table("tmp_peak_count_read1b",header=FALSE)
    peak_count_read2b = read.table("tmp_peak_count_read2b",header=FALSE)
    merge_common_only_count_read1a = read.table("tmp_merge_common_read1a", header=FALSE)
    merge_common_only_count_read1b = read.table("tmp_merge_common_read1b", header=FALSE)
    merge_common_only_count_read2a = read.table("tmp_merge_common_read2a", header=FALSE)
    merge_common_only_count_read2b = read.table("tmp_merge_common_read2b", header=FALSE)
    merge_common_peak_count_read1a = read.table("tmp_merge_common_peak_count_read1a",header=FALSE)
    merge_common_peak_count_read2a = read.table("tmp_merge_common_peak_count_read2a",header=FALSE)
    merge_common_peak_count_read1b = read.table("tmp_merge_common_peak_count_read1b",header=FALSE)
    merge_common_peak_count_read2b = read.table("tmp_merge_common_peak_count_read2b",header=FALSE)
    common_count_mat = cbind(common_peak_count_read1a[,4], common_peak_count_read1b[,4], 
        common_peak_count_read2a[,4], common_peak_count_read2b[,4])
    all_count_mat = cbind(peak_count_read1a[,4], peak_count_read1b[,4], 
        peak_count_read2a[,4], peak_count_read2b[,4])
    common_merge_count_mat = cbind(merge_common_only_count_read1a[,4], merge_common_only_count_read1b[,4],
        merge_common_only_count_read2a[,4], merge_common_only_count_read2b[,4])
    all_merge_count_mat = cbind(merge_common_peak_count_read1a[,4], merge_common_peak_count_read1b[,4], 
        merge_common_peak_count_read2a[,4], merge_common_peak_count_read2b[,4])
} else { 
    common_peak_count_read1 = read.table("tmp_common_peak_count_read1",header=FALSE)
    common_peak_count_read2 = read.table("tmp_common_peak_count_read2",header=FALSE)
    peak_count_read1 = read.table("tmp_peak_count_read1",header=FALSE)
    peak_count_read2 = read.table("tmp_peak_count_read2",header=FALSE)
    merge_common_only_count_read1 = read.table("tmp_merge_common_read1", header=FALSE)
    merge_common_only_count_read2 = read.table("tmp_merge_common_read2", header=FALSE)
    merge_common_peak_count_read1 = read.table("tmp_merge_common_peak_count_read1",header=FALSE)
    merge_common_peak_count_read2 = read.table("tmp_merge_common_peak_count_read2",header=FALSE)
    common_count_mat = cbind(common_peak_count_read1, common_peak_count_read2)
    all_count_mat = cbind(peak_count_read1, peak_count2)
    common_merge_count_mat = cbind(merge_common_only_count_read1, merge_common_only_count_read2)
    all_merge_count_mat = cbind(merge_common_peak_count_read1, merge_common_peak_count_read2)
}
table_MA = read.table("tmp_MAnorm.bed",header=FALSE)
table_merge_MA = read.table("tmp_MAnorm_merge.bed",header=FALSE)

#####################################################################################################
# new way + rescaling, works for 1 replicate only (0 replicates will give rowMeans error)
#M = log2(rowMeans(common_count_mat[,1:(ncol(common_count_mat)/2)]+1)/
#    rowMeans(common_count_mat[,((ncol(common_count_mat)/2)+1):ncol(common_count_mat)]+1))
#M = log2((rowMeans(common_count_mat[,1:2])+1)/(rowMeans(common_count_mat[,3:4])+1)) #same as above
#A = rowMeans(log2(common_count_mat+1)) #same as below
#see rest of plotMA code
# pretty always returns a 2+5 vector
#smoothScatter(A,M,cex=1,main="MA plot before rescaling (common peaks)",  xlim=pretty(A)[c(1, 7)])

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
normalized_all_count_mat = normalizeMA(common_count_mat+1, all_count_mat+1, method="rlm")-1
#subset_mat = common_count_mat+1; full_mat = all_count_mat+1 #DEBUG

png(paste0(f_prefix, '_MAsmooth_before_rescaling.png'))
plotSmoothMA(common_count_mat)
dev.off()

png(paste0(f_prefix, '_MAsmooth_before_rescaling_all.png'))
plotSmoothMA(all_count_mat)
dev.off()

png(paste0(f_prefix, '_MAsmooth_before_rescaling_all_common.png'))
plotSmoothMA(all_count_mat, sel = table_MA[,4] == "common_peak1" | table_MA[,4] == "common_peak2")
dev.off()



require(qvalue)
# DIFFBIND WAY (only works for replicates (group = 4))
require(edgeR)
res = DGEList(normalized_all_count_mat)
res$samples$group = c(1,1,2,2)
res$counts = round(res$counts) #method requires integers
res = calcNormFactors(res,method="TMM") #possibly not needed
res$design = model.matrix(~res$samples$group)
res = estimateGLMCommonDisp(res,res$design)
# when this is run without estimateGLMTrendedDisp called, it will 
#  "squeeze tagwise dispersions towards common dispersion" edgeR 2.8.2
# this means that when using plotBCV it wont show the trend line
res = estimateGLMTagwiseDisp(res,res$design)  
res$GLM = glmFit(res,res$design)
res$LRT = glmLRT(res$GLM,2)
#out = topTags(res$LRT n = nrow(normalized_all_count_mat)) #res$LRT$table stores all the fun stuff

png(paste0(f_prefix, '_MAsmooth_after_rescaling.png'))
#plotSmoothMA(normalized_all_count_mat)
plotSmoothMA(res$LRT)
dev.off()

res$LRT$qv = qvalue(res$LRT$table$PValue)
png(paste0(f_prefix, '_MAsmooth_after_rescaling_sig.png'))
#plotSmoothMA(normalized_all_count_mat, pvals=res$LRT$qv$qvalues)
plotSmoothMA(res$LRT, pvals=res$LRT$qv$qvalues)
dev.off()

# do the same thing for merged dataset
colnames(common_merge_count_mat) = c(x_name, y_name)
colnames(all_merge_count_mat) = c(x_name, y_name)
normalized_all_merge_count_mat = normalizeMA(common_merge_count_mat+1, all_merge_count_mat+1, method="rlm")-1

png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling.png'))
plotSmoothMA(common_merge_count_mat)
dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling_all.png'))
plotSmoothMA(all_merge_count_mat)
dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_before_rescaling_all_common.png'))
plotSmoothMA(all_merge_count_mat, sel = table_merge_MA[,4] == "merged_common_peak")
dev.off()


#DIFFBIND
resm = DGEList(normalized_all_merge_count_mat)
resm$samples$group = c(1,1,2,2)
resm$counts = round(resm$counts) #method requires integers
resm = calcNormFactors(resm,method="TMM") #possibly not needed
resm$design = model.matrix(~resm$samples$group)
resm = estimateGLMCommonDisp(resm,resm$design)
# estimateGLMCommonDisp() calculates logCPM abundance and common dispersion
# logCPM is around log2(rowMeans(cpm(resm))) (did not look up exact calculation)
resm = estimateGLMTagwiseDisp(resm,resm$design)  
resm$GLM = glmFit(resm,resm$design)
resm$LRT = glmLRT(resm$GLM,2)

resm$LRT$qv = qvalue(resm$LRT$table$PValue)

png(paste0(f_prefix, '_MAsmooth_merge_after_rescaling.png'))
plotSmoothMA(resm$LRT)
dev.off()

png(paste0(f_prefix, '_MAsmooth_merge_after_rescaling_sig.png'))
plotSmoothMA(resm$LRT, pvals=resm$LRT$qv$qvalues)
dev.off()


if(0) { # old way
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
table_MA[,13] = res$LRT$table$logFC
table_MA[,14] = res$LRT$table$logCPM
table_MA[,15] = -log10(res$LRT$table$PValue)
table_MA[,16] = -log10(res$LRT$qv$qvalues)

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
#careful about using qvalues package because diff pi0 can affect
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
table_merge_MA[,13] = resm$LRT$table$logFC
table_merge_MA[,14] = resm$LRT$table$logCPM
table_merge_MA[,15] = -log10(resm$LRT$table$PValue)
table_merge_MA[,16] = -log10(resm$LRT$qv$qvalues)

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

write.table(table_merge_MA,"MAnorm_result_commonPeak_merged.xls",sep="\t",quote=FALSE,row.names=FALSE)

rm(table_MA, table_merge_MA)
source("../edgeRmatricies.R")
