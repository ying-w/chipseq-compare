#slight modifications to MAnorm to be faster
#library(robustreg) #required for rlm
#####################################################################################################
# ideally we would just pass in a boolean subset index instead of full matrix
# but because of the way that the files are generated in earlier step, this is easier
# code mostly copied from affy::normalize.loess function (shown below)
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
                message("Done with",j,"vs",k,"in iteration",iter,"\n")
            } #for k
        } #for j
        full_mat = full_mat - means
        change = max(colMeans((means[subsample,])^2))

        if(verbose)
        message("Iteration: ", iter, ", change: ", change,"\n")
    } #while

    if ((change > epsilon) & (maxit > 1))
    warning(paste("No convergence after", maxit, "iterations.\n"))

    if(log.it) {
        return(2^full_mat)
    } else
    return(full_mat)
}
# modified from https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/affy/R/normalize.loess.R
# code is a bit cleaner than loess.normalize.R
normalize.loess = function(mat, subsample=sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
epsilon=10^-2, maxit=1, log.it=TRUE, verbose=TRUE, span=2/3,
family.loess="symmetric"){

    J = dim(mat)[2]
    II = dim(mat)[1]
    if(log.it){
        mat = log2(mat)
    }

    change = epsilon +1
    iter = 0
    w = c(0, rep(1,length(subsample)), 0) ##this way we give 0 weight to the
    ##extremes added so that we can interpolate

    while(iter < maxit){
        iter = iter + 1
        means = matrix(0,II,J) ##contains temp of what we substract

        for (j in 1:(J-1)){
            for (k in (j+1):J){
                y = mat[,j] - mat[,k]
                x = (mat[,j] + mat[,k]) / 2
                index = c(order(x)[1], subsample, order(-x)[1])
                ##put endpoints in so we can interpolate
                xx = x[index]
                yy = y[index]
                aux = loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
                aux = predict(aux, data.frame(xx=x)) / J
                means[, j] = means[, j] + aux
                means[, k] = means[, k] - aux
                if (verbose)
                cat("Done with",j,"vs",k,"in iteration",iter,"\n")
            }
        }
        mat = mat - means
        change = max(colMeans((means[subsample,])^2))

        if(verbose)
        cat(iter, change,"\n")
    }

    if ((change > epsilon) & (maxit > 1))
    warning(paste("No convergence after", maxit, "iterations.\n"))

    if(log.it) {
        return(2^mat)
    } else
    return(mat)
}

# function for calculating pvalue -> use with edgeR
# https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/DiffBind/R/analyze.R
# pv.DEedgeR function
# note: there is an 'edgeRGLM' option in DBA, not sure what it does though
# called from DBA.R as 
# dba.plotMA = function(DBA, contrast=1, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE, fold=0, bNormalized=TRUE,
#                     factor="", bXY=FALSE, dotSize=.33, bSignificant=TRUE, bSmooth=TRUE, ...)
#   res = pv.DBAplotMA(DBA, contrast=contrast, method=method, bMA=!bXY, bXY=bXY, th=th, bUsePval=bUsePval, fold=fold,
#                      facname=factor, bNormalized=bNormalized, cex=dotSize, 
#                      bSignificant = bSignificant, bSmooth=bSmooth,  ...)
pv.DBAplotMA = function(pv,contrast,method='edgeR',bMA=T,bXY=F,th=0.1,bUsePval=F,fold=0,facname="",bNormalized=T,
cex=.15,bSignificant=T, bSmooth=T,...) {

    if(missing(contrast)){
        contrast=1:length(pv$contrasts)
    } else {
        if(contrast > length(pv$contrasts)) {
            stop('Specified contrast number is greater than number of contrasts')
            return(NULL)
        }
    }

    plotfun = plot
    if (bSmooth) {
        plotfun = smoothScatter
    }


    numSites = nrow(pv$vectors)

    for(con in 1:length(contrast)) {
        conrec = pv$contrasts[[contrast[con]]]
        for(meth in method) {
            res = pv.DBAreport(pv,contrast=contrast[con],method=meth,bUsePval=T,th=100,bNormalized=bNormalized)
            if(!is.null(res)) {
                if(bUsePval) {
                    idx = res$"p-value" <= th
                    tstr = "p"
                } else {
                    idx = res$FDR <= th
                    tstr = "FDR"
                }
                idx = idx & (abs(res$Fold) >= fold)
                if(bMA){
                    xmin  = floor(min(res$Conc))
                    xmax  = ceiling(max(res$Conc))
                    ymin  = floor(min(res$Fold))
                    ymax  = ceiling(max(res$Fold))
                    if(bSmooth | !bSignificant) {
                        ## !!!! THIS IS THE KEY FUNCTION !!!!
                        ## cex = 0.33, facname = "" (factor name), tstr = FDR
                        plotfun(res$Conc,res$Fold,pch=20,cex=cex,
                            xaxp=c(xmin,xmax,xmax-xmin),xlim=c(xmin,xmax),
                            xlab='log concentration',
                            yaxp=c(ymin,ymax,(ymax-ymin)),ylim=c(ymin,ymax),
                            ylab=sprintf('log fold change: %s - %s',conrec$name1,conrec$name2),
                            main=sprintf('%s Binding Affinity: %s vs. %s (%s %s < %1.3f)',
                            facname, conrec$name1,conrec$name2,sum(idx),tstr,th),...)                  
                    } else {
                        plotfun(res$Conc[!idx],res$Fold[!idx],pch=20,cex=cex,
                        xaxp=c(xmin,xmax,xmax-xmin),xlim=c(xmin,xmax),
                        xlab='log concentration',
                        yaxp=c(ymin,ymax,(ymax-ymin)),ylim=c(ymin,ymax),
                        ylab=sprintf('log fold change: %s - %s',conrec$name1,conrec$name2),
                        main=sprintf('%s Binding Affinity: %s vs. %s (%s %s < %1.3f)',
                        facname, conrec$name1,conrec$name2,sum(idx),tstr,th),...)
                    }
                    if(bSignificant) { ## COLOR SIGNIFICANT POINTS
                        points(res$Conc[idx],res$Fold[idx],pch=20,cex=cex,col=2)
                    }
                    abline(h=0,col='dodgerblue') ## DRAW HORIZONTAL
                } # if bMA
                if(bXY){
                    xmin  = floor(min(res[,5]))
                    xmax  = ceiling(max(res[,5]))
                    ymin  = floor(min(res[,6]))
                    ymax  = ceiling(max(res[,6]))
                    xymin = min(xmin,ymin)
                    xymin = max(xymin,0)
                    xymax = max(xmax,ymax)
                    plotfun(res[!idx,6],res[!idx,5],pch=20,cex=cex,col=1,
                    xaxp=c(xymin,xymax,xymax-xymin),xlim=c(xymin,xymax),
                    xlab=sprintf('log concentration :%s',conrec$name2),
                    yaxp=c(xymin,xymax,(xymax-xymin)),ylim=c(xymin,xymax),
                    ylab=sprintf('log concentration :%s',conrec$name1),
                    main=sprintf('%s Binding Affinity: %s vs. %s (%s %s < %1.3f)',
                    facname, conrec$name1,conrec$name2,sum(idx),tstr,th),...)
                    points(res[idx,6],res[idx,5],pch=20,cex=cex,col=2)
                    abline(0,1,col='dodgerblue')
                } #if bXY
            } #if(!is.null(res)
        } 
    }    
}
# mergeReps can be append or add, c("append", "add") gives comparison warning
plotSmoothMA = function(mat, plotfun = smoothScatter, mergeReps = "append", sel=rep(FALSE, nrow(mat))) {
    #takes UNlogged mat
    #TODO replace X & Y w/colnames
    if(ncol(mat) > 2) {
        if(mergeReps == "append") { 
            message("Appending replicates")
            x = log2(mat[,c(1:(ncol(mat)/2))]+1)
            y = log2(mat[,c((ncol(mat)/2+1):ncol(mat))]+1)
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

    # M = res$Fold
    # A = res$Conc
    # idx = res$FDR <= th
    plotfun(A, M, pch = 20, cex = 0.33, xlim = c(0, ceiling(max(A))),
        xlab='log concentration',
        ylab=sprintf('log fold change: %s - %s',"x","y"),
        main=sprintf('%s Binding Affinity: %s vs. %s (%s of %s highlighted w/%s < %1.3f)',
            'TF', 'x', 'y', sum(sel), length(A), "pval",1)) #if length(A) is null then nothing will be printed
    abline(h=0,col='dodgerblue')
    points(A[sel], M[sel], pch=20, cex=0.33,, col=2) #highlight significant
}

# when called by dba.analyze()
#  bSubControl = TRUE, bFullLibrarySize = FALSE, bTagwise = TRUE, blockList is NULL
# when using debug(DiffBind:::pv.DEedgeR) remmeber to set DBA$config$RunParallel = FALSE or debugger will hang
pv.DEedgeR = function(pv,group1,group2,label1="Group 1",label2="Group 2",blockList=NULL,
                bSubControl=F,bFullLibrarySize=F,bTagwise=T,bGLM=T,bNormOnly=F) {

    fdebug('Enter pv.DEedgeR')

    #require(edgeR)
    res = pv.DEinit(pv,group1,group2,label1,label2,method='edgeR',
    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize)
    res = calcNormFactors(res,method="TMM")
    fdebug(sprintf('calcNormFactors: %f',res$counts[7,1]))

    if(bNormOnly) {
        return(res)    
    }

    if(is.null(blockList)) {
        fdebug('blockList is NULL')
    } else {
        fdebug('blockList is not NULL')
    }
    fdebug('pv.DEedgeR: check for blocking factor')
    if(is.null(blockList)) {
        fdebug('pv.DEedgeR: NO blocking factor')
        res = estimateCommonDisp(res)
        fdebug(sprintf('estimateCommonDisp: %f',res$counts[7,1]))
        if(bGLM){
            res$design = model.matrix(~res$samples$group)
            if(bTagwise) {
                res = estimateGLMCommonDisp(res,res$design)
                res = estimateGLMTagwiseDisp(res,res$design)
            } else {
                res = estimateGLMCommonDisp(res,res$design)
            }
            res$GLM = glmFit(res,res$design)
            res$LRT = glmLRT(res$GLM,2)
        } else {
            if(bTagwise){
                res = estimateTagwiseDisp(res,prior.df=50,trend="none")
                #res = estimateTagwiseDisp(res,prior.n=getPriorN(res),trend="movingave")
                res$db     = exactTest(res,dispersion='tagwise')
            } else {
                res$db     = exactTest(res,dispersion='common')    
            }
        }
        fdebug(sprintf('Fit and test: %f',res$counts[7,1]))
        fdebug('pv.DEedgeR: estimateTagwiseDisp complete')
        
        fdebug(sprintf('pv.DEedgeR: exactTest complete:%s-%s',res$db$comparison[1],res$db$comparison[2]))
        #res$db$fdr = topTags(res$db,nrow(res$db$counts))
    } else {
        fdebug('pv.DEedgeR: BLOCKING FACTOR')
        
        targets = pv.blockFactors(pv,group1,group2,label1,label2,blockList)
        if(is.null(targets)){
            return(res)    
        }
        res$samples = data.frame(cbind(targets,res$samples[,2:3]))
        
        attr =  blockList[[1]]$attribute
        if(attr=='Replicate') {
            res$designmatrix = model.matrix(~ Replicate + group,data = targets)
        } else if(attr=='Tissue') {
            res$designmatrix = model.matrix(~ Tissue + group,data = targets)
        } else if(attr=='Factor') {
            res$designmatrix = model.matrix(~ Factor + group,data = targets)
        } else if(attr=='Condition') {
            res$designmatrix = model.matrix(~ Condition + group,data = targets)
        } else if(attr=='Caller') {
            res$designmatrix = model.matrix(~ Caller + group,data = targets)
        } else if(attr=='Treatment') {
            res$designmatrix = model.matrix(~ Treatment + group,data = targets)
        } else if(attr=='Block') {
            res$designmatrix = model.matrix(~ Block + group,data = targets)
        } else {
            warning('Unsupported blocking attribute: ',attr,call.=FALSE)
            return(NULL)    
        }
        message('edgeR multi-factor analysis.')
        res = calcNormFactors(res)
        res = estimateGLMCommonDisp(res,res$designmatrix)
        if(bTagwise) {
            res = estimateGLMTagwiseDisp(res,res$designmatrix)
        }
        res$GLM = glmFit(res,res$designmatrix)
        res$LRT = glmLRT(res$GLM,ncol(res$designmatrix))
        res$counts=NULL     
        #res$fdr = topTags(res$LRT,nrow(res$counts))
    }

    res$bSubControl      = bSubControl
    res$bFullLibrarySize = bFullLibrarySize

    fdebug(sprintf('Exit pv.DEedgeR: %f',res$counts[7,1]))
    return(res)    

}
# context for above call
# for(i in 1:length(pv$contrast)) {     
    # res = pv.DEedgeR(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
        # pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
        # bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,bGLM=bGLM)

    # if(!is.null(pv$contrasts[[i]]$blocklist)) {
            # res$block = pv.DEedgeR(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
            # pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
            # pv$contrasts[[i]]$blocklist,
            # bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise)   
    # }
    # reslist = pv.listadd(reslist,res)   
# }
#####################################################################################################
## no more functions, lets begin
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
    merge_common_peak_count_read1a = read.table("tmp_merge_common_peak_count_read1a",header=FALSE)
    merge_common_peak_count_read2a = read.table("tmp_merge_common_peak_count_read2a",header=FALSE)
    merge_common_peak_count_read1b = read.table("tmp_merge_common_peak_count_read1b",header=FALSE)
    merge_common_peak_count_read2b = read.table("tmp_merge_common_peak_count_read2b",header=FALSE)
    common_count_mat = cbind(common_peak_count_read1a[,4], common_peak_count_read1b[,4], 
        common_peak_count_read2a[,4], common_peak_count_read2b[,4])
    all_count_mat = cbind(peak_count_read1a[,4], peak_count_read1b[,4], 
        peak_count_read2a[,4], peak_count_read2b[,4])
} else { 
    common_peak_count_read1 = read.table("tmp_common_peak_count_read1",header=FALSE)
    common_peak_count_read2 = read.table("tmp_common_peak_count_read2",header=FALSE)
    peak_count_read1 = read.table("tmp_peak_count_read1",header=FALSE)
    peak_count_read2 = read.table("tmp_peak_count_read2",header=FALSE)
    merge_common_peak_count_read1 = read.table("tmp_merge_common_peak_count_read1",header=FALSE)
    merge_common_peak_count_read2 = read.table("tmp_merge_common_peak_count_read2",header=FALSE)
    common_count_mat = cbind(common_peak_count_read1, common_peak_count_read2)
    all_count_mat = cbind(peak_count_read1, peak_count2)
}
table_MA = read.table("tmp_MAnorm.bed",header=FALSE)
table_merge_MA = read.table("tmp_MAnorm_merge.bed",header=FALSE)

#####################################################################################################
# new way + rescaling, works for 1 replicate only (0 replicates will give rowMeans error)
#M = log2(rowMeans(common_count_mat[,1:(ncol(common_count_mat)/2)]+1)/
#    rowMeans(common_count_mat[,((ncol(common_count_mat)/2)+1):ncol(common_count_mat)]+1))
M = log2((rowMeans(common_count_mat[,1:2])+1)/(rowMeans(common_count_mat[,3:4])+1)) #same as above
A = rowMeans(log2(common_count_mat+1)) #same as below
#see rest of plotMA code
# pretty always returns a 2+5 vector
#smoothScatter(A,M,cex=1,main="MA plot before rescaling (common peaks)",  xlim=pretty(A)[c(1, 7)])

normalized_all_count_mat = normalizeMA(common_count_mat+1, all_count_mat+1, method="rlm")
#subset_mat = common_count_mat+1
#full_mat = all_count_mat+1

png('MAsmooth_before_rescaling_all.png')
plotSmoothMA(all_count_mat)
dev.off()

png('MAsmooth_before_rescaling.png')
plotSmoothMA(common_count_mat)
dev.off()

png('MAsmooth_after_rescaling.png')
plotSmoothMA(normalized_all_count_mat)
dev.off()

#run edgeR code below
png('MAsmooth_after_rescaling_sig.png')
plotSmoothMA(normalized_all_count_mat, sel=res$LRT$table$PValue<0.01)
dev.off()

# old way
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


#edgeR example
# group <- factor(c(1,1,2,2))
# y <- DGEList(counts=x,group=group)
# y <- calcNormFactors(y)
# y <- estimateCommonDisp(y)
# y <- estimateTagwiseDisp(y)
# et <- exactTest(y)
# topTags(et)

# design <- model.matrix(~group)
# y <- estimateGLMCommonDisp(y,design)
# y <- estimateGLMTrendedDisp(y,design)
# y <- estimateGLMTagwiseDisp(y,design)
# fit <- glmFit(y,design)
# lrt <- glmLRT(fit,coef=2)
# topTags(lrt)

# DIFFBIND WAY
require(edgeR)
res = DGEList(normalized_all_count_mat)
res$samples$group = c(1,1,2,2)
res$counts = round(res$counts) #remove warnings
res = calcNormFactors(res,method="TMM")
res$design = model.matrix(~res$samples$group)
res = estimateGLMCommonDisp(res,res$design)
# when this is run without estimateGLMTrendedDisp called, it will 
#  "squeeze tagwise dispersions towards common dispersion" edgeR 2.8.2
res = estimateGLMTagwiseDisp(res,res$design)  
res$GLM = glmFit(res,res$design)
res$LRT = glmLRT(res$GLM,2)
out = topTags(res$LRT n = nrow(normalized_all_count_mat)) #res$LRT$table stores all the fun stuff
# eventually gets stored in GR$contrasts[[1]]$edgeR$



table_MA[,5] = peak_count_read1
table_MA[,6] = peak_count_read2
table_MA[,7] = M_rescaled
table_MA[,8] = A_rescaled
table_MA =as.data.frame(table_MA)
table_MA[,9] = 0
log2_peak_count_read1_rescaled = as.matrix(log2_peak_count_read1_rescaled)
peak_count_read2 = as.matrix(peak_count_read2)
for (n in c(1:nrow(table_MA))) {
#        cat(n,'\t',round(2^log2_peak_count_read1_rescaled[n]),'\t',peak_count_read2[n],'\n')
    table_MA[n,9]=-log10(pval(round(2^(log2_peak_count_read1_rescaled[n])),peak_count_read2[n]))
}


colnames(table_MA)[1] = "chr"
colnames(table_MA)[2] = "start"
colnames(table_MA)[3] = "end"
colnames(table_MA)[4] = "description"
colnames(table_MA)[5] = "#raw_read_1"
colnames(table_MA)[6] = "#raw_read_2"
colnames(table_MA)[7] = "M_value_rescaled"
colnames(table_MA)[8] = "A_value_rescaled"
colnames(table_MA)[9] = "-log10(p-value)"

write.table(table_MA,"MAnorm_result.xls",sep="\t",quote=FALSE,row.names=FALSE)

# table_merge
table_merge_MA[,5] = merge_common_peak_count_read1
table_merge_MA[,6] = merge_common_peak_count_read2
table_merge_MA[,7] = merge_M_rescaled
table_merge_MA[,8] = merge_A_rescaled
table_merge_MA =as.data.frame(table_merge_MA)
table_merge_MA[,9] = 0
log2_merge_common_peak_count_read1_rescaled = as.matrix(log2_merge_common_peak_count_read1_rescaled)
merge_common_peak_count_read2 = as.matrix(merge_common_peak_count_read2)
for (n in c(1:nrow(table_merge_MA))) {
#        cat(n,'\t',round(2^log2_merge_common_peak_count_read1_rescaled[n]),'\t',merge_common_peak_count_read2[n],'\n')
    table_merge_MA[n,9]=-log10(pval(round(2^(log2_merge_common_peak_count_read1_rescaled[n])),merge_common_peak_count_read2[n]))
}


colnames(table_merge_MA)[1] = "chr"
colnames(table_merge_MA)[2] = "start"
colnames(table_merge_MA)[3] = "end"
colnames(table_merge_MA)[4] = "description"
colnames(table_merge_MA)[5] = "#raw_read_1"
colnames(table_merge_MA)[6] = "#raw_read_2"
colnames(table_merge_MA)[7] = "M_value_rescaled"
colnames(table_merge_MA)[8] = "A_value_rescaled"
colnames(table_merge_MA)[9] = "-log10(p-value)"

write.table(table_merge_MA,"MAnorm_result_commonPeak_merged.xls",sep="\t",quote=FALSE,row.names=FALSE)
