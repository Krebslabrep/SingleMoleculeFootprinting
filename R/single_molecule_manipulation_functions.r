#' Generic function to design bins based on different user case
#' 
#' @param RegionsOfInterest GRanges of length == 1 (for single TFs and promoters) or >= 1 for TF clusters containing the biological loci to sort around
#' @param BinType one of "singleTF", "TFcluster", "PromoterDM", "PromoterMM". Defaults to "singleTF"
#' @param userBin custom user-passed bins in the format list(bin1 = c(start, end), bin2 = c(start, end), ...). Defaults to NULL. If used the arg BinType is ignored
#' @param StrandAware TRUE or FALSE (default). If BinType is "singleTF" or "TFcluster", StrandAware argument is ignored
#' 
#' @import GenomicRanges
#' 
#' @export
#' 
MakeBins <- function(RegionsOfInterest, BinType = "singleTF", userBin = NULL, StrandAware = FALSE){
    
  ### Select type of bin and bin coordinates
  ### Select viewpoint (start coordinate of the region of interest)
  message("Designing sorting bins")
  
  if(!is.null(userBin)){
    # user-defined bin
    bins <- userBin
  } else if(BinType == "singleTF"){
    viewPoint = start(RegionsOfInterest) + (end(RegionsOfInterest)-start(RegionsOfInterest))/2
    bins <- list(up = c(-35,-25), TFBS = c(-15,15), down = c(25,35))
  } else if(BinType == "TFcluster"){
    RegionsOfInterest = sort(RegionsOfInterest, by = ~ seqnames + start + end)
    viewPoint = start(RegionsOfInterest) + (end(RegionsOfInterest)-start(RegionsOfInterest))/2
    bins = list(up = c(-35,-25), TFBS = c(-7,7), down = c(25,35))
  } else if(BinType == "PromoterDM"){
    viewPoint <- as.numeric(start(RegionsOfInterest))
    bins <- list(up = c(-58, -43), TATA = c(-36, -22), INR = c(-6, 14), DPE = c(28, 47))
  } else if(BinType == "PromoterMM"){
    viewPoint <- as.numeric(start(RegionsOfInterest))
    bins <- list(up = c(-58, -43), TATA = c(-36, -22), PolII = c(28, 47), down=c(54, 69))
  }
  
  if(BinType == "singleTF"){
    
    binList = GRanges(seqnames = seqnames(RegionsOfInterest),
                      ranges = IRanges(
                        start = c(viewPoint+bins[[1]][1], viewPoint+bins[[2]][1], viewPoint+bins[[3]][1]),
                        end = c(viewPoint+bins[[1]][2], viewPoint+bins[[2]][2], viewPoint+bins[[3]][2])),
                      strand = strand(RegionsOfInterest),
                      seqinfo = seqinfo(RegionsOfInterest))
   
  } else if(BinType == "TFcluster"){
 
    binList = IRanges(start = c(min(viewPoint)+bins[[1]][1], viewPoint+bins[[2]][1], max(viewPoint)+bins[[3]][1]),
                      end = c(min(viewPoint)+bins[[1]][2], viewPoint+bins[[2]][2], max(viewPoint)+bins[[3]][2]))
    
  } else {
    
    ### bin coordinates
    binList <- lapply(1:length(bins), function(i){
      # select bin
      binAbsPos <- bins[[i]]
      binName <- names(bins[i])
      # get bin according to StrandAware status
      if(StrandAware == TRUE){
        BinsCoordinates <- GRanges(seqnames = seqnames(RegionsOfInterest),
                                   ranges = IRanges(
                                     start = ifelse(as.logical(strand(RegionsOfInterest) == '+'), viewPoint + binAbsPos[1], viewPoint - binAbsPos[2]),
                                     end = ifelse(as.logical(strand(RegionsOfInterest) == '+'), viewPoint + binAbsPos[2], viewPoint - binAbsPos[1])),
                                   strand = strand(RegionsOfInterest),
                                   seqinfo = seqinfo(RegionsOfInterest))
      } else if(StrandAware == FALSE){
        BinsCoordinates <- GRanges(seqnames = seqnames(RegionsOfInterest),
                                   ranges = IRanges(
                                     start =  viewPoint + binAbsPos[1],
                                     end = viewPoint + binAbsPos[2]),
                                   seqinfo = seqinfo(RegionsOfInterest))
      }
      # bin name (eg. upstream, TFBS, etc.)
      BinsCoordinates$name <- binName
      # feature name
      names(BinsCoordinates) <- names(RegionsOfInterest)
      return(BinsCoordinates)
  })
    binList <- unique(sort(unlist(GRangesList(binList))))
  }
  
  return(binList)
}

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
#' binMethylationValues = BinMethylation(MethSM = MethSM, Bin = BinsCoordinates[1])
#'
BinMethylation = function(MethSM, Bin){

  binCytosines = colnames(MethSM)[as.numeric(colnames(MethSM)) >= start(Bin) & as.numeric(colnames(MethSM)) <= end(Bin)]

  # Summarise methylation status of each read
  if (length(binCytosines) >= 1){
    binSummarisedMeth = round(rowMeans_drop0(MethSM[,binCytosines,drop=FALSE]) - 1)
    binSummarisedMeth = binSummarisedMeth[!(is.na(binSummarisedMeth))]
    return(binSummarisedMeth)
  } else if (length(binCytosines) == 0){
    message(paste0("!!!     [", start(Bin), ";", end(Bin), "]", " bin overlaps with no covered Cytosines   !!!"))
    return(NA)
  }



}

#' Sort reads by single TF
#'
#' @param MethSM Single molecule matrix
#' @param BinsCoordinates IRanges object of absolute coordinates for sorting bins
#' @param coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed
#' @param strand add strand awareness to single molecule sorting. Defaults "+"
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
SortReads = function(MethSM, BinsCoordinates, coverage=NULL, strand = "+"){

  message("Collecting summarized methylation for bins")
  binMethylationList = lapply(seq_along(BinsCoordinates), function(i){
    BinMethylation(MethSM = MethSM, Bin = BinsCoordinates[i])
  })

	message("Subsetting those reads that cover all bins")
	ReadsSubset = Reduce(intersect, lapply(binMethylationList, function(x){names(x)}))
	if (length(ReadsSubset) < coverage){
	  message(paste0("Less than ", coverage, " reads found to cover all sorting bins...skipping"))
	  return(list())
	}

	message("Summarizing reads into patterns")
	binMethylationList_subset = lapply(binMethylationList, function(x){as.character(x[ReadsSubset])})
	MethPattern = Reduce(paste0, binMethylationList_subset)
	
	# reverse pattern if strand == "-"
	if(strand == "-"){
	  MethPattern <- unlist(lapply(MethPattern, stri_reverse))
	}

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
#' @param bins [DEPRECATED] list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-15,15), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the TFBS.
#'             bins[[2]] represents the TFBS bin, with coordinates relative to the center of the TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the TFBS.
#' @param coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed. Defaults to 20
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
SortReadsBySingleTF = function(MethSM, TFBS, bins = list(c(-35,-25), c(-15,15), c(25,35)), coverage = 20){

  message("Designing sorting bins")
  warning("Argument bins DEPRECATED. Now hardcoded in MakeBins function")
  # TFBS_center = start(TFBS) + (end(TFBS)-start(TFBS))/2
  # BinsCoordinates = IRanges(start = c(TFBS_center+bins[[1]][1], TFBS_center+bins[[2]][1], TFBS_center+bins[[3]][1]),
  #                           end = c(TFBS_center+bins[[1]][2], TFBS_center+bins[[2]][2], TFBS_center+bins[[3]][2]))
  BinsCoordinates = MakeBins(RegionsOfInterest = TFBS, BinType = "singleTF", userBin = NULL, StrandAware = FALSE)
  SortedReads = lapply(MethSM, SortReads, BinsCoordinates = BinsCoordinates, coverage = coverage)
  
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
#' @param coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed. Defaults to 30
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
SortReadsByTFCluster = function(MethSM, TFBS_cluster, bins = list(c(-35,-25), c(-7,7), c(25,35)), coverage = 30){

  # message("Sorting TFBSs by genomic coordinates")
  # TFBS_cluster = sort(TFBS_cluster, by = ~ seqnames + start + end)
  message("Designing sorting bins")
  warning("Argument bins DEPRECATED. Now hardcoded in MakeBins function")
  # TFBS_centers = start(TFBS_cluster) + (end(TFBS_cluster)-start(TFBS_cluster))/2
  # BinsCoordinates = IRanges(start = c(min(TFBS_centers)+bins[[1]][1], TFBS_centers+bins[[2]][1], max(TFBS_centers)+bins[[3]][1]),
  #                     end = c(min(TFBS_centers)+bins[[1]][2], TFBS_centers+bins[[2]][2], max(TFBS_centers)+bins[[3]][2]))
  BinsCoordinates = MakeBins(RegionsOfInterest = TFBS_cluster, BinType = "TFcluster", userBin = NULL, StrandAware = FALSE)
  SortedReads = lapply(MethSM, SortReads, BinsCoordinates = BinsCoordinates, coverage = coverage)
  
  return(SortedReads)

}

#' Wrapper to SortReads for Promoter case
#'
#' @param MethSM Single molecule matrix list as returned by CallContextMethylation
#' @param TSSs Transcription Start Sites to use for sorting, passed as a GRanges object of length >= 1
#' @param bins list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-7,7), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the most upstream TFBS.
#'             bins[[2]] represents all the TFBS bins, with coordinates relative to the center of each TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the most downstream TFBS.
#' @param coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed. Defaults to 30
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
SortReadsByPromoter = function(MethSM, TFBS_cluster, bins = list(c(-35,-25), c(-7,7), c(25,35)), coverage = 30){
  
  message("Sorting TFBSs by genomic coordinates")
  TFBS_cluster = sort(TFBS_cluster, by = ~ seqnames + start + end)
  message("Designing sorting bins")
  TFBS_centers = start(TFBS_cluster) + (end(TFBS_cluster)-start(TFBS_cluster))/2
  BinsCoordinates = IRanges(start = c(min(TFBS_centers)+bins[[1]][1], TFBS_centers+bins[[2]][1], max(TFBS_centers)+bins[[3]][1]),
                            end = c(min(TFBS_centers)+bins[[1]][2], TFBS_centers+bins[[2]][2], max(TFBS_centers)+bins[[3]][2]))
  
  SortedReads = lapply(MethSM, SortReads, BinsCoordinates = BinsCoordinates, coverage = coverage)
  
  return(SortedReads)
  
}

#' Convenience for calculating state frequencies
#' 
#' @param SortedReads List of sorted reads (can be multiple samples) as returned by either read sorting function (SortReads, SortReadsBySingleTF, SortReadsByTFCluster)
#' @param states states reporting the biological interpretation of patterns as return by either OneTFstates or TFpairStates functions. If NULL (default) will return frequencies without biological interpretation.
#' 
#' @importFrom tibble tibble
#' 
#' @return tibble with state frequency information
#' 
#' @export
#' 
StateQuantification = function(SortedReads, states){
  
  if(all(isEmpty(SortedReads))){
    
    return(tibble(Sample = NA, State = NA, Counts = NA, Freqs = NA))
    
  }
  
  if (is.null(states)){
    
    Patterns = unique(unlist(lapply(SortedReads, names)))
    states = split(Patterns, Patterns)
    
  }
  
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
