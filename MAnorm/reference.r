# some R functions from bioconductor svn that I used as reference
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

#In diffbind data is eventually stored in GR$contrasts[[1]]$edgeR$
