#slight modifications to MAnorm to be faster
#library(robustreg) #required for rlm
library(MASS) #required for rlm
#####################################################################################################
# modified from https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/affy/R/normalize.loess.R
# code is a bit cleaner than loess.normalize.R
normalize.loess <- function(mat, subset=sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
                            epsilon=10^-2, maxit=1, log.it=TRUE, verbose=TRUE, span=2/3,
                            family.loess="symmetric"){

  J <- dim(mat)[2]
  II <- dim(mat)[1]
  if(log.it){
    mat <- log2(mat)
  }

  change <- epsilon +1
  iter <- 0
  w <- c(0, rep(1,length(subset)), 0) ##this way we give 0 weight to the
                                      ##extremes added so that we can interpolate

  while(iter < maxit){
    iter <- iter + 1
    means <- matrix(0,II,J) ##contains temp of what we substract

    for (j in 1:(J-1)){
      for (k in (j+1):J){
        y <- mat[,j] - mat[,k]
        x <- (mat[,j] + mat[,k]) / 2
        index <- c(order(x)[1], subset, order(-x)[1])
        ##put endpoints in so we can interpolate
        xx <- x[index]
        yy <- y[index]
        aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
        aux <- predict(aux, data.frame(xx=x)) / J
        means[, j] <- means[, j] + aux
        means[, k] <- means[, k] - aux
        if (verbose)
          cat("Done with",j,"vs",k,"in iteration",iter,"\n")
      }
    }
    mat <- mat - means
    change <- max(colMeans((means[subset,])^2))

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

# modified from https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/DiffBind/R/analyze.R
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
              if(bSignificant) {
                 points(res$Conc[idx],res$Fold[idx],pch=20,cex=cex,col=2)
              }
              abline(h=0,col='dodgerblue')
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
#####################################################################################################

if(length(list.files(pattern="read1a.bed")) == 1) {  #11 arguments used
	
} else { 
	common_peak_count_read1 <- read.table("tmp_common_peak_count_read1",header=FALSE)
	common_peak_count_read2 <- read.table("tmp_common_peak_count_read2",header=FALSE)
	peak_count_read1 <- read.table("tmp_peak_count_read1",header=FALSE)
	peak_count_read2 <- read.table("tmp_peak_count_read2",header=FALSE)
	merge_common_peak_count_read1 <- read.table("tmp_merge_common_peak_count_read1",header=FALSE)
	merge_common_peak_count_read2 <- read.table("tmp_merge_common_peak_count_read2",header=FALSE)
}
table_MA <-read.table("tmp_MAnorm.bed",header=FALSE)
table_merge_MA <- read.table("tmp_MAnorm_merge.bed",header=FALSE)

M <- as.matrix(log2((common_peak_count_read1+1)/(common_peak_count_read2+1)))
A <- as.matrix(0.5*log2((common_peak_count_read1+1)*(common_peak_count_read2+1)))

linear <- lm(M~A)$coefficients
#b<-lm(M~A)$coefficients
#b<-robustRegBS(M,A,beta=linear)
b<-rlm(M~A)$coefficients

png('MAplot_before_rescaling.png')
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
ma.plot(A,M,cex=1,main="MA plot before rescaling (common peaks)")
abline(b,col="green")
dev.off()

cat("M = b[1] + b[2] * A\n")
log2_peak_count_read1 <- log2(peak_count_read1 + 1)
log2_peak_count_read2 <- log2(peak_count_read2 + 1)
log2_peak_count_read1_rescaled <- (2-b[2])*log2_peak_count_read1/(2+b[2]) - 2*b[1]/(2+b[2]);
#you get the above forumla after solving for read1 when you plug in M and A
M_rescaled <- (log2_peak_count_read1_rescaled - log2_peak_count_read2);
A_rescaled <- (log2_peak_count_read1_rescaled + log2_peak_count_read2)/2;

png('MAplot_after_rescaling.png')
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
ma.plot(as.matrix(A_rescaled),as.matrix(M_rescaled),cex=1,main=" MA plot after rescaling (all peaks)")
dev.off ()

log2_merge_common_peak_count_read1 <- log2(merge_common_peak_count_read1 + 1)
log2_merge_common_peak_count_read2 <- log2(merge_common_peak_count_read2 + 1)
log2_merge_common_peak_count_read1_rescaled <- (2-b[2])*log2_merge_common_peak_count_read1/(2+b[2]) - 2*b[1]/(2+b[2]);
merge_M_rescaled <- (log2_merge_common_peak_count_read1_rescaled - log2_merge_common_peak_count_read2);
merge_A_rescaled <- (log2_merge_common_peak_count_read1_rescaled + log2_merge_common_peak_count_read2)/2;

# function for calculating pvalue -> use with edgeR

table_MA[,5] <- peak_count_read1
table_MA[,6] <- peak_count_read2
table_MA[,7] <- M_rescaled
table_MA[,8] <- A_rescaled
table_MA <-as.data.frame(table_MA)
table_MA[,9] <- 0
log2_peak_count_read1_rescaled <- as.matrix(log2_peak_count_read1_rescaled)
peak_count_read2 <- as.matrix(peak_count_read2)
for (n in c(1:nrow(table_MA))) {
#        cat(n,'\t',round(2^log2_peak_count_read1_rescaled[n]),'\t',peak_count_read2[n],'\n')
        table_MA[n,9]<--log10(pval(round(2^(log2_peak_count_read1_rescaled[n])),peak_count_read2[n]))
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
table_merge_MA[,5] <- merge_common_peak_count_read1
table_merge_MA[,6] <- merge_common_peak_count_read2
table_merge_MA[,7] <- merge_M_rescaled
table_merge_MA[,8] <- merge_A_rescaled
table_merge_MA <-as.data.frame(table_merge_MA)
table_merge_MA[,9] <- 0
log2_merge_common_peak_count_read1_rescaled <- as.matrix(log2_merge_common_peak_count_read1_rescaled)
merge_common_peak_count_read2 <- as.matrix(merge_common_peak_count_read2)
for (n in c(1:nrow(table_merge_MA))) {
#        cat(n,'\t',round(2^log2_merge_common_peak_count_read1_rescaled[n]),'\t',merge_common_peak_count_read2[n],'\n')
        table_merge_MA[n,9]<--log10(pval(round(2^(log2_merge_common_peak_count_read1_rescaled[n])),merge_common_peak_count_read2[n]))
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
