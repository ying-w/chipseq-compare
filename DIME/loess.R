#DIME_preprocess.R
winsize = 500
samp = "high1"
refc = "none1"
chrlist = c(1:22,"M","X","Y")

for(i in chrlist){
	print(paste("Processing chr", i, sep=""))

	file_sample = paste("./data/", samp, "/", paste(samp, winsize, paste("chr",i,sep=""),sep="_"),sep="")
	file_reference = paste("./data/", refc, "/", paste(refc, winsize, paste("chr",i,sep=""),sep="_"),sep="")
	
	sami = read.table(file_sample)
	refi = read.table(file_reference)
	
	diffCt = sami[,1]-refi[,1]
	meanCt = (sami[,1]+refi[,1])/2
	chr_size = nrow(refi)
	
	spanMean = 0.6
	spanVar = 0.1
	
	pc = proc.time()
	#Warning: NaNs
	diffCt_loesstmp = loess(diffCt ~ meanCt, span=spanMean ) #takes a long time
	#diffCt_loessMean = predict(diffCt_loesstmp, meanCt[,1]) 
	diffCt_loessMean = predict(diffCt_loesstmp)
	diffCt_loessMean[is.na(diffCt_loessMean)] = 0
	norm_diffCt = diffCt - diffCt_loessMean
	norm_diffCt_loesstmp = loess(abs(norm_diffCt) ~ meanCt, span=spanVar)
	#norm_diffCt_loessVar = predict(norm_diffCt_loesstmp, meanCt)
	norm_diffCt_loessVar = predict(norm_diffCt_loesstmp)
	
	#eps = 2^(-52)
	#norm_diffCt[norm_diffCt==0] = eps;
	#norm_diffCt_loessVar[norm_diffCt_loessVar==0] = eps;
	
	varNorm_diffCt = norm_diffCt/norm_diffCt_loessVar
	varNorm_diffCt[is.na(varNorm_diffCt)] = 0
	
	outfile = paste("./data/", paste("vs", samp, refc, winsize, paste("chr",i,sep=""),sep="_"),sep="")
	write.table(varNorm_diffCt, file=outfile, sep="\t", quote=F, row.names=F, col.names=F)
	print(proc.time()-pc)
	#library(DIME); hn = read.table(_); hn.d = DIME(hn[,1]) #takes a while 
}

#
