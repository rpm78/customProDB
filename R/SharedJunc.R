##' Load multiple BED files and output a GRange object with junctions present in multiple samples. 
##'
##' This function allows to limit junctions that are present in at least m out of n BED files.
##' @title Generate shared junctions dataset from multiple BED files
##' @param juns a list of GRanges object which input from multiple VCF files using function InputVcf.
##' @param share_num Junctions must occurs in this number of samples to be consider. Two options, percentage format or sample number.
##' @param ... additional arguments
##' @return a GRange object that contains the shared junctions
##' @author Xiaojing Wang
##' @examples
##' 
##' path <- system.file("extdata/beds", package="customProDB")
##' bedFiles<- paste(path, '/', list.files(path, pattern="*bed$"), sep='')
##' juncs <- lapply(bedFiles, function(x) Bed2Range(x, skip=1, covfilter=5))
##' shared <- SharedJunc(juncs, share_num=2)
##' shared
##'

SharedJunc <-  function(juns, share_num=2, 
	#ext_up=100, ext_down=100, 
	...)
    {
   
        jungls <- GRangesList(juns)
        
        jungls_basic <- lapply(jungls, function(x){
                c(x, ignore.mcols=TRUE)
                })
        jungls_basic <- GRangesList(jungls_basic)
        total <- unlist(jungls_basic)       
        
		if(grepl('%', share_num)){
            share_num <- round(as.numeric(gsub('%', '', share_num)) * 
                                length(juns) / 100) 
        }else share_num <- as.numeric(share_num)
        
		uniquetotal <- unique(total)

        ctab <- do.call('cbind', lapply(jungls, function(x) 
                            countOverlaps(uniquetotal, x, type='equal')))
			
        sumcount <- apply(ctab, 1, sum)
        index <- which(sumcount >= share_num)
        sharejunc <- uniquetotal[index]
        
		#### index matrix of sharejunc in jungls
		tt <- lapply(jungls, function(x) 
				data.frame(findOverlaps(sharejunc, x, type='equal')))
		sharejuncIDXjungls <- matrix(NA, length(sharejunc), length(jungls))
		for(x in 1:length(tt)){
			sharejuncIDXjungls[tt[[x]][, 1], x] <- tt[[x]][, 2]
		}
		
		###max part1 length
		part1_jungls <-  do.call('cbind', lapply(1:length(jungls), function(x) 
				data.frame(jungls[[x]])[sharejuncIDXjungls[, x], 'part1_len']))
		part1_jungls[is.na(part1_jungls)] <- 0
		part1_len <- apply(part1_jungls, 1, max)	
		
		###max part2 length
		part2_jungls <-  do.call('cbind', lapply(1:length(jungls), function(x) 
				data.frame(jungls[[x]])[sharejuncIDXjungls[, x], 'part2_len']))
		part2_jungls[is.na(part2_jungls)] <- 0
		part2_len <- apply(part2_jungls, 1, max)	
		
		###mean coverage across samples have that junction (round)
		cov_jungls <-  do.call('cbind', lapply(1:length(jungls), function(x) 
				data.frame(jungls[[x]])[sharejuncIDXjungls[, x], 'cov']))
		cov_mean <- round(apply(cov_jungls, 1,  function(x) mean(x, na.rm=TRUE)))
		
		
        junRange <- GRanges(seqnames=seqnames(sharejunc), ranges=ranges(sharejunc), 
                strand=strand(sharejunc), id=paste('JUNC', 1:length(index), sep=''), 
                cov=cov_mean, part1_len=part1_len, part2_len=part2_len, 
                part1_sta=start(sharejunc)-part1_len+1, part1_end=start(sharejunc), 
                part2_sta=end(sharejunc), part2_end=end(sharejunc)+part2_len+1)
 
}