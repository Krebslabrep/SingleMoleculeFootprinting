#' Summarize methylation inside sorting bins
#'
#' @param MethSM Single molecule matrix
#' @param TFBS Transcription factor binding site to use for sorting, passed as a GRanges object of length 1
#' @param bin vector of two integers representing the coordinate of a bin relative to the center of the TFBS
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#'
#' @return Reads covering bin with their summarized methylation status
#'
#' @export
#'
#' @examples
#'
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#'
#' binMethylationValues = BinMethylation(MethSM = MethSM, TFBS = TFBSs, bin = c(-15,15))
#'
BinMethylation = function(MethSM, TFBS, bin){

  midP = start(resize(TFBS, 1, fix='center'))
  binRange = GRanges(seqnames(TFBS), IRanges(start = midP+bin[1], end = midP+bin[2]))

  # Find Cs of interest
  binCytosines = colnames(MethSM)[as.numeric(colnames(MethSM)) >= start(binRange) & as.numeric(colnames(MethSM)) <= end(binRange)]

  # Summarise methylation status of each read
  if (length(binCytosines) > 1){
    binSummarisedMeth = round(rowMeans(MethSM[,binCytosines], na.rm = TRUE) + 0.01) # adding tipping over value of 0.01 so that corner cases where true mean == 0.5 will round to 1
  } else if (length(binCytosines) == 1){
    binSummarisedMeth = MethSM[,binCytosines]
  } else if (length(binCytosines) == 0){
    stop(paste0("[", bin[1], ";", bin[2], "]", " bin overlaps with no covered Cytosines"))
  }

  binSummarisedMeth = binSummarisedMeth[!(is.na(binSummarisedMeth))]

  return(binSummarisedMeth)

}

#' Sort reads by single TF
#'
#' @param MethSM Single molecule matrix
#' @param TFBS Transcription factor binding site to use for sorting
#' @param BinsCoord list of 3 bin coordinates relative to the center of the TFBS.
#' @param SortByCluster T/F
#'
#' @import BiocGenerics
#'
#' @return list of sorted reads
#'
#' @export
#'
#' @examples
#'
#' BinsCoord = list(c(-35,-25), c(-15,15), c(25,35))
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#'
#' SortedReads = SortReads(MethSM, TFBSs, BinsCoord, SortByCluster = FALSE)
#'
SortReads = function(MethSM, TFBS, BinsCoord, SortByCluster){

  message("Collecting summarized methylation for bins")
  if(!SortByCluster){
    message("Single TF mode")
    binMethylationList = lapply(BinsCoord, BinMethylation, MethSM = MethSM, TFBS = TFBS)
  } else {
    # The following is terrible code and should be fixed
    message("TF cluster mode")
    TFBSs = sort(TFBS)
    UpstreamBinMethylation = BinMethylation(MethSM, TFBSs[1], BinsCoord[[1]])
    TFBSsBinMethylation = lapply(seq_along(TFBSs), function(i){BinMethylation(MethSM, TFBSs[i], BinsCoord[[2]])}) # apply to each TFBS
    DownstramBinMethylation = BinMethylation(MethSM, TFBSs[length(TFBSs)], BinsCoord[[3]])
    binMethylationList = list()
    binMethylationList[[1]] = UpstreamBinMethylation
    for (i in seq_along(TFBSsBinMethylation)){
      binMethylationList[[1+i]] = TFBSsBinMethylation[[i]]
    }
    binMethylationList[[length(binMethylationList)+1]] = DownstramBinMethylation
  }

	# Subset the reads that cover all bins
	ReadsSubset = Reduce(intersect, lapply(binMethylationList, function(x){names(x)}))

	# Get pattern for each read
	binMethylationList_subset = lapply(binMethylationList, function(x){as.character(x[ReadsSubset])})
	MethPattern = Reduce(paste0, binMethylationList_subset)

	# Split reads by pattern
	if(length(ReadsSubset)>0){
	  sortedReadslist = split(ReadsSubset, MethPattern) #paste(names(TFBS), MethPattern, sep='_'))
	}else{
	  sortedReadslist = list()
	  }

	return(sortedReadslist)

}

#' Wrapper to SortReads for single TF case
#'
#' @param MethSM Single molecule matrix
#' @param TFBS Transcription factor binding site to use for sorting, passed as a GRanges object of length 1
#'
#' @return List of reads sorted by single TF
#'
#' @export
#'
#' @examples
#'
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#'
#' SortedReads = SortReadsBySingleTF(MethSM = MethSM, TFBS = TFBSs)
#'
SortReadsBySingleTF = function(MethSM, TFBS){

  BinsCoord = list(c(-35,-25), c(-15,15), c(25,35))
  SortedReads = SortReads(MethSM, TFBS, BinsCoord, SortByCluster = FALSE)
  return(SortedReads)

}

#' Wrapper to SortReads for TF cluster case
#'
#' @param MethSM Single molecule matrix
#' @param TFBSs Transcription factor binding sites to use for sorting, passed as a GRanges object of length > 1
#'
#' @return List of reads sorted by TF cluster
#'
#' @export
#'
#' @examples
#'
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106216, 88106253), c(88106226, 88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1", "NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305215, 4305216)))
#'
#' SortedReads = SortReadsByTFCluster(MethSM = MethSM, TFBSs = TFBSs)
#'
SortReadsByTFCluster = function(MethSM, TFBSs){

  BinsCoord = list(c(-35,-25), c(-7,7), c(25,35))
  SortedReads = SortReads(MethSM, TFBS = TFBSs, BinsCoord, SortByCluster = TRUE)
  return(SortedReads)

}
