#' Function to extract regions of interest from a tabix-indexed bed.gz file
#' 
#' @param file Path the bgzipped bed file to parse
#' @param regions GRanges or string (chr2L:2000-5000) with the coordinates of the regions to extract
#' @param strict restrict to bed lines within the regions of interest
#' @param col_names string of the column names contained in the bed file (automatic naming of the 3 first column)
#' 
#' @import GenomicRanges
#' @import tibble
#' @importFrom Rsamtools scanTabix
#' @importFrom dplyr distinct bind_rows
#' @importFrom readr type_convert
#' 
#' @return Return a tibble or dataframe containing the bed entries mapped only for the regions of interest
#'
#' @export
#'
#' @examples
#' Bedfile = system.file("extdata", "methBedDemo.bed.gz", package = "SingleMoleculeFootprinting", mustWork = T)
#' range =  GRanges(c("chr2R", "chr3L"), IRanges(c(6250000, 10000000), width=10000000))
#' mbedColnames = c("chrom","start","end","readname","mstring","scores","context")
#' region.tb <- getFromTabix(file=Bedfile, regions=range, col_names=mbedColnames, tibble=FALSE)
#'
getFromTabix <- function(file, regions, strict=FALSE, col_names=NULL, tibble=FALSE){

	###### NEED TO BE REMOVED AT DEPLOYMENT ############
	dep.l <- c('GenomicRanges', 'tibble', 'dplyr', 'Rsamtools', 'readr')
	for (p in dep.l) {if(!suppressMessages(require(p, character.only = TRUE))) stop(paste0("Package not found: ", p))}
	###### NEED TO BE REMOVED AT DEPLOYMENT ############

	# extract region of interest using tabix
	errMsg <- "The regions argument syntax isn't correct!\nMust be like: 'chr2L:1-3000,chr3R:400-8000'\n"
	if (!class(regions) == "GRanges"){
		if(length(grep(":", regions)) == 0 | length(grep("-", regions)) == 0) {stop(errMsg)}
		splitregion <- unlist(strsplit(regions, ','))
		regions <- unlist(GRangesList(lapply(splitregion, function(r){
			tmp = unlist(strsplit(r, ":"))
			coord = as.numeric(unlist(strsplit(tmp[2], "-")))
			return(GRanges(tmp[1], IRanges(start=coord[1], end=coord[2]), strand='*'))
		})))
	}
	reads.list = Rsamtools::scanTabix(file, param=regions)
	ncolBED <- length(unlist(strsplit(reads.list[[1]][1], "\t")))
	reads.df <- do.call(dplyr::bind_rows, lapply(seq_along(reads.list), function(l){
		# conversion to tibble
			df <- tibble::as_tibble(do.call(rbind, strsplit(reads.list[[l]], "\t")), .name_repair = ~ paste0('X', seq(1:ncolBED)))
			df <- suppressMessages(readr::type_convert(df)) %>% dplyr::distinct()
			
		# filter for overlapping reads
			if(strict & dim(df)[1] > 0){
				df <- df %>% filter(V2 > start(regions)[l])  %>% filter(V3 < end(regions)[l])
			}
		return(df)
	}))
	# correct column names if necessary
	if(!is.null(col_names)){
		if(dim(reads.df)[1] == 0){
			reads.df <- as_tibble(matrix(nrow = 0, ncol = length(col_names)), .name_repair = ~ col_names)
		} else {
			colnames(reads.df)  <- col_names
		}
	} else {
		colnames(reads.df) <- c('chr', 'start', 'end', paste0('unknown_', seq(dim(reads.df)[2]-3)))
	}
	if(tibble){
		return(reads.df)
	} else {
		return(as.data.frame(reads.df))
	}
}


#' Function to get QuasR::qMeth() like output from single molecule methylation call from Nanopore data
#' 
#' BLABLAVLA
#' 
#' @param SM_files data frame ....
#' @param regExp 
#' @param flankedRegion genomic distance in base-pair to expand the regExp region single molecule methylation call. Defaults to 0.
#' @param mbedColnames string containing the headers of the single-molecule bed file.
#' @param modeGCG mode to keep ambiguous context methylation value. Will keep only the CpG (mode='cpg'), GpC (mode='gpc') or none (mode='strict') of methylation value at GCG context. Default to 'cpg'.
#' 
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges reduce
#' 
#' @return list with 4 elements:
#' - aid: readID
#' - Cid: Cytosine coordinates
#' - strand: context strand
#' - meth: cytosine methylation call (methylated=1, unmethylated=0, undetermined=NA)
#'
#' @examples
#'
#' Bedfile = system.file("extdata", "methBedDemo.bed.gz", package = "SingleMoleculeFootprinting", mustWork = T)
#' SM_files = data.frame(bed=filePath[grepl('_calls.tsv.bgz$', filePath)], tabix=filePath[grepl('_calls.tsv.bgz.tbi$', filePath)], row.names=names(pattern[k]))
#' QuasRprj = GetQuasRprj(Qinput, BSgenome.Mmusculus.UCSC.mm10)
#' Samples = QuasR::alignments(QuasRprj)[[1]]$SampleName
#' sample = Samples[1]
#' MethGR = QuasR::qMeth(QuasRprj[grep(sample, Samples)], mode="allC", range, collapseBySample = TRUE, keepZero = TRUE)
#'

.getQuasRqMeth <- function(SM_files, regExp, flankedRegion=0, mbedColnames=c("chrom","start","end","readname","mstring","scores","context"), modeGCG='cpg'){
	QuasRmat <- lapply(1:nrow(SM_files), function(i){
		region.tb <- getFromTabix(file=SM_files[i,1], regions=regExp, col_names=mbedColnames, tibble=FALSE, strict=FALSE)
		#region.tb <- region.tb[!duplicated(region.tb$readname),] ## supplementary alignments: same read could be find multiple times (removal random in regard of primary or secondary alignment)
		if(dim(region.tb)[1] == 0){
			return(list(data.frame(), GRanges()))
		} else {
			return(list(.mBedToQuasRmat(mBed=region.tb), GRanges(region.tb[,1:4])))
		}
		})
	readOrder <- sort(unlist(GenomicRanges::GRangesList(lapply(QuasRmat, function(i){i[[2]]}))))
	QuasRmat <- do.call(rbind, lapply(QuasRmat, function(i){i[[1]]}))
	if(nrow(QuasRmat) == 0){
		return(list(aid=c(), Cid=c(), strand=c(), meth=c()))
	} else {
		# order the matrix by reads mapping coord
			QuasRmat$aid <- factor(QuasRmat$aid, unique(readOrder$readname))
			QuasRmat <- QuasRmat[order(QuasRmat$aid, QuasRmat$Cid),]
			rownames(QuasRmat) <- seq(nrow(QuasRmat))
		# GCG are duplicated in case of separated context call
			QuasRmat <- .correctMatForGCGdup(mat=QuasRmat, mode=modeGCG)
		# reduce for only the windows of observation
			toKeep=QuasRmat$Cid >=  start(regExp) - flankedRegion & QuasRmat$Cid <=  end(regExp) + flankedRegion
			QuasRmat <- QuasRmat[toKeep,]
		return(list(aid=as.character(QuasRmat$aid),
				Cid=QuasRmat$Cid,
				strand=QuasRmat$strand,
				meth=QuasRmat$meth))
	}
}

#' Function that converts single-molecule methylation string (contained in Nanopolish bed output files) to single-molecule cytosine methylation data frame
#' 
#' @param mBed data frame containing methylBed reads out of Nanopore methylation call
#'
#' @return data frame with 4 columns:
#' - ReadID (aid)
#' - Cytosine coordinate (Cid)
#' - strand (strand)
#' - meth: translated single molecule methylation string into binary methylation (0=unmethylated, 1=methylated, NA=undetermined)
#'
#' @examples
#'
#' Bedfile = system.file("extdata", "methBedDemo.bed.gz", package = "SingleMoleculeFootprinting", mustWork = T)
#' range =  GRanges(c("chr2R", "chr3L"), IRanges(c(6250000, 10000000), width=10000000))
#' mbedColnames = c("chrom","start","end","readname","mstring","scores","context")
#' region.tb <- getFromTabix(file=Bedfile, regions=range, col_names=mbedColnames, tibble=FALSE)
#' QuasRmat <- .mBedToQuasRmat(mBed=region.tb)
#'
.mBedToQuasRmat <- function(mBed){
	mat.l <- lapply(1:nrow(mBed), function(r){
		read <- mBed[r,]
		mstring =read$mstring
		call = strsplit(mstring, "[0-9]+")[[1]][-1]
		pos = cumsum(as.numeric(strsplit(mstring, "[a-z]")[[1]]))
		mat <- data.frame(
			aid=rep(read$readname, length(pos)),
			Cid=as.numeric(read$start) + pos + 1,
			strand=rep('+', length(pos)),
			meth=NA)
		mat$meth[which(call == "m")] = 1
		mat$meth[which(call == "u")] = 0
		return(mat)
	})
	return(do.call(rbind, mat.l))
}

#' Function to fix GCG context duplicates in single-molecule cytosine methylation data frame
#'
#' If multiple call found for one cytosine on the same molecule (not possible with our improved Nanopolish version):
#' 	First, if one of them has undetermined status correct with the other methylation call.
#' 	Second, keep only the CpG or GpC methylation value (mode='cpg' or 'gpc' respectively).
#' 	If the mode is 'strict', all ambiguous cytosine in GCG with 2 methylation call are removed
#' 
#' @param mat single-molecule cytosine methylation data frame
#' @param mode mode to keep ambiguous context methylation value {cpg|gpc|strict}
#'
#' @return GCG corrected single-molecule cytosine methylation data frame with 4 columns:
#' - ReadID (aid)
#' - Cytosine coordinate (Cid)
#' - strand (strand)
#' - meth: translated single molecule methylation string into binary methylation (0=unmethylated, 1=methylated, NA=undetermined)
#'
#' @examples
#'
#' Bedfile = system.file("extdata", "methBedDemo.bed.gz", package = "SingleMoleculeFootprinting", mustWork = T)
#' range =  GRanges(c("chr2R", "chr3L"), IRanges(c(6250000, 10000000), width=10000000))
#' mbedColnames = c("chrom","start","end","readname","mstring","scores","context")
#' region.tb <- getFromTabix(file=Bedfile, regions=range, col_names=mbedColnames, tibble=FALSE)
#' QuasRmat <- .mBedToQuasRmat(mBed=region.tb)
#' QuasRmat.GCGcorreted <- .correctMatForGCGdup(mat=QuasRmat, mode='cpg')
#'
.correctMatForGCGdup <- function(mat, mode='cpg'){
	dup <- as.numeric(rownames(mat[duplicated(paste(mat$aid, mat$Cid, sep='_')),])) #should be produced by the gpc call
	if(length(dup) > 0){
		# correct NAs with the other call if possible
		summary(is.na(mat[dup-1, 'meth']) == is.na(mat[dup, 'meth']))
		mat[dup-1, 'meth'] <- ifelse(is.na(mat[dup-1, 'meth']), mat[dup, 'meth'], mat[dup-1, 'meth'])
		mat[dup, 'meth'] <- ifelse(is.na(mat[dup, 'meth']), mat[dup-1, 'meth'], mat[dup, 'meth'])

		# correct value that are differential called between the cpg and gpc
		mat[dup-1, 'meth'] <- ifelse(
			Vectorize(isTRUE)(mat[dup, 'meth'] == mat[dup-1, 'meth']),
			mat[dup, 'meth'], 
			if(mode == 'cpg'){
				mat[dup-1, 'meth']
			} else if(mode == 'gpc'){
				mat[dup, 'meth']
			} else if(mode == 'strict'){
				NA
			} else {
				stop('The mode of .correctMatForGCG() is not correct!\n')
			})
		return(mat[-dup,])
	} else {
		return(mat)
	}
}

