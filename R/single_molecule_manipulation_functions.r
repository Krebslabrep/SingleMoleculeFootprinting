#' Summarize methylation inside sorting bins
#'
#' @param MethSM Single molecule matrix
#' @param Bin IRanges object with absolute coordinates for single sorting bin.
#'
#' @import GenomicRanges
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
#' TFBS_center = start(TFBS) + (end(TFBS)-start(TFBS))/2
#' BinsCoordinates = IRanges(start = c(TFBS_center+bins[[1]][1], TFBS_center+bins[[2]][1], TFBS_center+bins[[3]][1]),
#'                           end = c(TFBS_center+bins[[1]][2], TFBS_center+bins[[2]][2], TFBS_center+bins[[3]][2]))
#'
#' binMethylationValues = BinMethylation(MethSM = MethSM, Bin = BinsCoordinates[1]))
#'
BinMethylation = function(MethSM, Bin){

  binCytosines = colnames(MethSM)[as.numeric(colnames(MethSM)) >= start(Bin) & as.numeric(colnames(MethSM)) <= end(Bin)]

  # Summarise methylation status of each read
  if (length(binCytosines) > 1){
    binSummarisedMeth = round(rowMeans(MethSM[,binCytosines], na.rm = TRUE))
  } else if (length(binCytosines) == 1){
    binSummarisedMeth = MethSM[,binCytosines[1]]
  } else if (length(binCytosines) == 0){
    stop(paste0("[", start(Bin), ";", end(Bin), "]", " bin overlaps with no covered Cytosines"))
  }

  binSummarisedMeth = binSummarisedMeth[!(is.na(binSummarisedMeth))]

  return(binSummarisedMeth)

}

#' Sort reads by single TF
#'
#' @param MethSM Single molecule matrix
#' @param BinsCoordinates IRanges object of absolute coordinates for sorting bins
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
#' TFBS_center = start(TFBS) + (end(TFBS)-start(TFBS))/2
#' BinsCoordinates = IRanges(start = c(TFBS_center+bins[[1]][1], TFBS_center+bins[[2]][1], TFBS_center+bins[[3]][1]),
#'                           end = c(TFBS_center+bins[[1]][2], TFBS_center+bins[[2]][2], TFBS_center+bins[[3]][2]))
#'
#' SortedReads = SortReads(MethSM, BinsCoordinates)
#'
SortReads = function(MethSM, BinsCoordinates){

  message("Collecting summarized methylation for bins")
  binMethylationList = lapply(seq_along(BinsCoordinates), function(i){
    BinMethylation(MethSM = MethSM, Bin = BinsCoordinates[i])
  })

	message("Subsetting those reads that cover all bins")
	ReadsSubset = Reduce(intersect, lapply(binMethylationList, function(x){names(x)}))

	message("Summarizing reads into patterns")
	binMethylationList_subset = lapply(binMethylationList, function(x){as.character(x[ReadsSubset])})
	MethPattern = Reduce(paste0, binMethylationList_subset)

	message("Splitting reads by pattern")
	if(length(ReadsSubset)>0){
	  sortedReadslist = split(ReadsSubset, MethPattern) #paste(names(TFBS), MethPattern, sep='_'))
	}else{
	  sortedReadslist = list()
	  }

	return(sortedReadslist)

}

#' Wrapper to SortReads for single TF case
#'
#' @param MethSM Single molecule matrix list as returned by CallContextMethylation
#' @param TFBS Transcription factor binding site to use for sorting, passed as a GRanges object of length 1
#' @param bins list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-15,15), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the TFBS.
#'             bins[[2]] represents the TFBS bin, with coordinates relative to the center of the TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the TFBS.
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
SortReadsBySingleTF = function(MethSM, TFBS, bins = list(c(-35,-25), c(-15,15), c(25,35))){

  message("Designing sorting bins")
  TFBS_center = start(TFBS) + (end(TFBS)-start(TFBS))/2
  BinsCoordinates = IRanges(start = c(TFBS_center+bins[[1]][1], TFBS_center+bins[[2]][1], TFBS_center+bins[[3]][1]),
                            end = c(TFBS_center+bins[[1]][2], TFBS_center+bins[[2]][2], TFBS_center+bins[[3]][2]))

  SortedReads = lapply(MethSM, SortReads, BinsCoordinates = BinsCoordinates)
  
  return(SortedReads)

}

#' Wrapper to SortReads for TF cluster case
#'
#' @param MethSM Single molecule matrix list as returned by CallContextMethylation
#' @param TFBS_cluster Transcription factor binding sites to use for sorting, passed as a GRanges object of length > 1
#' @param bins list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-7,7), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the most upstream TFBS.
#'             bins[[2]] represents all the TFBS bins, with coordinates relative to the center of each TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the most downstream TFBS.
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
#' SortedReads = SortReadsByTFCluster(MethSM = MethSM, TFBSs = TFBS_cluster)
#'
SortReadsByTFCluster = function(MethSM, TFBS_cluster, bins = list(c(-35,-25), c(-7,7), c(25,35))){

  message("Sorting TFBSs by genomic coordinates")
  TFBS_cluster = sort(TFBS_cluster, by = ~ seqnames + start + end)
  message("Designing sorting bins")
  TFBS_centers = start(TFBS_cluster) + (end(TFBS_cluster)-start(TFBS_cluster))/2
  BinsCoordinates = IRanges(start = c(min(TFBS_centers)+bins[[1]][1], TFBS_centers+bins[[2]][1], max(TFBS_centers)+bins[[3]][1]),
                      end = c(min(TFBS_centers)+bins[[1]][2], TFBS_centers+bins[[2]][2], max(TFBS_centers)+bins[[3]][2]))

  SortedReads = lapply(MethSM, SortReads, BinsCoordinates = BinsCoordinates)
  
  return(SortedReads)

}

#' Convenience for calculating state frequencies
#' 
#' @param SortedReads List of sorted reads (can be multiple samples) as returned by either read sorting function (SortReads, SortReadsBySingleTF, SortReadsByTFCluster)
#' @param states states reporting the biological interpretation of patterns as return by either OneTFstates or TFpairStates functions.
#' 
#' @return tibble with state frequency information
#' 
#' @export
#' 
StateQuantification = function(SortedReads, states){
  
  OrderedReads = lapply(SortedReads, function(sR){sR[as.character(unlist(states))]})
  
  Reduce(rbind,
  lapply(seq_along(OrderedReads), function(i){
    
    unlist(lapply(seq_along(states), function(j){
      length(unlist(OrderedReads[[i]][states[[j]]]))
    })) -> Counts
    
    tibble(Sample = names(OrderedReads)[i], 
           State = names(states), 
           Counts = Counts, 
           Freqs = (Counts/sum(Counts))*100)
    
  })) -> StateQuantification_tbl
  
  return(StateQuantification_tbl)
  
}

#' Convenience for calculating state frequencies after sorting reads by single TF
#' 
#' wraps around StateQuantification function
#' 
#' @param SortedReads List of sorted reads (can be multiple samples) as returned by SortReadsBySingleTF (or SortReads run with analogous parameters)
#' 
#' @return tibble with state frequency information
#' 
#' @export
#' 
StateQuantificationBySingleTF = function(SortedReads){
  
  states = OneTFstates()
  res = StateQuantification(SortedReads = SortedReads, states = states)
  return(res)
  
}

#' Convenience for calculating state frequencies after sorting reads by TF pair
#' 
#' wraps around StateQuantification function
#' 
#' @param SortedReads List of sorted reads (can be multiple samples) as returned by SortReadsByTFCluster run for clusters of size 2 (or SortReads run with analogous parameters)
#' 
#' @return tibble with state frequency information
#' 
#' @export
#' 
StateQuantificationByTFPair = function(SortedReads){
  
  states = TFpairStates()
  res = StateQuantification(SortedReads = SortedReads, states = states)
  return(res)
  
}
