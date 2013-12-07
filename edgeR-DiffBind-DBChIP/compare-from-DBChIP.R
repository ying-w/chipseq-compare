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