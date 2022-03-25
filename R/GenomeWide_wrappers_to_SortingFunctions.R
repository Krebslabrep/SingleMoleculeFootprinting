#' Convenience function to arrange a list of given TFBSs into clusters
#' 
#' For each TFBS, the genomic neighborhood defined by max_cluster_width will be scanned for adjacent TFBSs. 
#' The hits will be filtered for min_intersite_distance where, in case of overlapping TFBSs, the second TFBS will be arbitrarily dropped.
#' These TFBSs plus the central "anchoring" one will define a TFBS cluster. 
#' This approach implies that the same TFBS can be employed to design multiple clusters in a sliding-window fashion.
#' 
#' @param TFBSs GRanges object of TFBSs
#' @param max_intersite_distance maximum allowed distance in base pairs between two TFBS centers for them to be considered part of the same cluster. Defaults to 75.
#' @param min_intersite_distance minimum allowed distance in base pairs between two TFBS centers for them not to be discarded as overlapping. 
#'                               This parameter should be set according to the width of the bins used for later sorting. Defaults to 15.
#' @param max_cluster_size maximum number of TFBSs to be contained in any given cluster. Defaults to 6
#' 
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges reduce
#' 
#' @return list with two elements: ClusterCoordinates (GRanges object of clusters coordinates) and ClusterComposition (GRangesList of sites for each cluster)
#' 
#' @export
#'
Arrange_TFBSs_clusters = function(TFBSs, max_intersite_distance = 75, min_intersite_distance = 15, max_cluster_size = 6){
  
  collection_window_width = max_intersite_distance*2
  TFBSs_resized_1 = resize(TFBSs,1,fix='center')
  TFBSs_resized_window = resize(TFBSs,collection_window_width,fix='center')
  
  Overlaps = findOverlaps(TFBSs_resized_window, TFBSs_resized_1, ignore.strand=TRUE)
  pair_dist = start(TFBSs_resized_1[subjectHits(Overlaps)]) - start(TFBSs_resized_1[queryHits(Overlaps)])
  
  message("Removing self-overlaps and redundant pairs")
  message("Removing pairs containing overlapping factors")
  Overlaps_filtered = Overlaps[pair_dist>min_intersite_distance,]
  
  message("Constructing GRanges object of clusters coordinates")
  TF_cluster = IRanges::reduce(GRanges(
    seqnames(TFBSs[queryHits(Overlaps_filtered)]),
    IRanges(start(TFBSs[queryHits(Overlaps_filtered)]), 
            end(TFBSs[subjectHits(Overlaps_filtered)]))
  )) # very rarely these exceed width of 300bp
  TF_cluster = sort(TF_cluster)
  
  message("Computing number of sites per cluster")
  TF_cluster$number_of_TF = countOverlaps(TF_cluster,TFBSs)
  message("Creating TFBS_cluster ID on the fly")
  names(TF_cluster) = paste0('TFBS_cluster_',seq_along(TF_cluster))
  
  message(paste0("Discaring clusters with more than ", max_cluster_size, "sites"))
  TF_cluster = TF_cluster[TF_cluster$numer_of_TF <= max_cluster_size]
  
  message("Constructing GRangesList of sites per cluster")
  Overlaps_clusters = findOverlaps(TF_cluster,TFBSs_resized_1)
  TF_list = split(TFBSs[subjectHits(Overlaps_clusters)], queryHits(Overlaps_clusters))
  names(TF_list) = names(TF_cluster)
  
  return(list(ClusterCoordinates = TF_cluster, ClusterComposition = TF_list))
  
}

#' Create methylation calling windows to call context methylation in one run for clusters lying proximally to each other
#' 
#' Relevant for genome-wide analyses
#' 
#' @param TFBS_cluster_coordinates TFBS cluster coordinates analogous to ClusterCoordinates object returned by Arrange_TFBSs_clusters function
#' @param max_intercluster_distance maximum distance between two consecutive TFBS clusters for them to be grouped in the same window
#' @param max_window_width upper limit to window width. This value should be adjusted according to the user's system as it determines the amount of memory used in the later context methylation call
#' @param min_cluster_width lower limit to window width. Corresponds to the scenario when a window contains a single TFBS cluster.
#' 
#' @import GenomicRanges
#' 
#' @return GRanges object of window coordinates to be used for more efficient calls of CallContextMethylation 
#' 
#' 
Create_MethylationCallingWindows = function(TFBS_cluster_coordinates, 
                                            max_intercluster_distance = 100000, 
                                            max_window_width = 5000000, 
                                            min_cluster_width = 600){
  
  message(paste0("Group TFBS_clusters that fall within ", max_intercluster_distance, "bp from each other in broader searching windows"))
  SearchingWindows = plyranges::reduce_ranges(resize(TFBS_cluster_coordinates, width = max_intercluster_distance, fix = 'center'))
  message("Trimming searching windows")
  start(SearchingWindows) = start(SearchingWindows) + ((max_intercluster_distance/2) - (min_cluster_width/2))
  end(SearchingWindows) = end(SearchingWindows) - ((max_intercluster_distance/2) - (min_cluster_width/2))
  
  TooLarge = width(SearchingWindows) > max_window_width
  
  while(sum(TooLarge) > 0){
    
    SmallerWindow = max_intercluster_distance*0.9
    message("Reducing too large searching windows")
    
    Overlaps = findOverlaps(SearchingWindows[TooLarge], TFBS_cluster_coordinates)
    SmallerWindows = plyranges::reduce_ranges(resize(TFBS_cluster_coordinates[subjectHits(Overlaps)], width = SmallerWindow, fix = "center"))
    SearchingWindows = sort(c(SearchingWindows[!TooLarge], SmallerWindows))
    
    max_intercluster_distance = SmallerWindow
    TooLarge = width(SearchingWindows) > max_window_width
    
  }
  
  return(sort(SearchingWindows))
  
}

#' Convenience wrapper to sort single molecule according to TFBS clusters at multiple sites in the genome
#' 
#' The function starts from a list of single TFBSs, arranges them into clusters, calls methylation at the interested sites and outputs sorted reads
#' 
#' @param sampleSheet QuasR pointer file
#' @param sample samples to use, from the SampleName field of the sampleSheet
#' @param genome BSgenome
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to 0.8. To skip this filtering step, set to NULL. For more information, check out the details section.
# #' @param clObj cluster object for parallel processing of multiple samples/RegionsOfInterest. For now only used by qMeth call for bulk methylation. Should be the output of a parallel::makeCluster() call
#' @param CytosinesToMask THIS IS DUCK-TAPE
#' @param TFBSs GRanges object of transcription factor binding sites coordinates
#' @param max_interTF_distance maximum distance between two consecutive TFBSs for them to be grouped in the same window
#' @param max_window_width upper limit to window width. This value should be adjusted according to the user's system as it determines the amount of memory used in the later context methylation call
#' @param min_cluster_width lower limit to window width. Corresponds to the scenario when a window contains a single TFBS.
#' @param bins list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-15,15), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the most upstream TFBS.
#'             bins[[2]] represents all the TFBS bins, with coordinates relative to the center of each TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the most downstream TFBS.
#' @param sorting_coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed. Defaults to 30.
#' @param cores number of cores to use for parallel processing of multiple Methylation Calling Windows (i.e. groupings of adjecent TFBS clusters)
#' 
#' @importFrom parallel mclapply
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' 
#' @return list where [[1]] is the TFBSs GRanges object describing coordinates TFBSs used to sort single molecules
#'                    [[2]] is a list of SortedReads nested per TFBS_cluster and sample
#'                    [[3]] is a tibble reporting the count (and frequency) of reads per state, sample and TFBS cluster
#' 
#' @export
#' 
SortReadsBySingleTF_MultiSiteWrapper = function(sampleSheet, sample, genome, coverage = 20, ConvRate.thr = 0.8, # clObj=NULL, ---> parameters passed to CallContextMethylation
                                                 CytosinesToMask = NULL,
                                                 TFBSs,
                                                 max_interTF_distance = 100000, max_window_width = 5000000, min_cluster_width = 600, # ---> parameters passed to Create_MethylationCallingWindows
                                                 sorting_coverage = 30, bins = list(c(-35,-25), c(-15,15), c(25,35)), # ---> parameters passed to SortReadsByTFCluster
                                                 cores = 1
){
  
  names(TFBSs) = paste0("TFBS_", seq(TFBSs))
  
  message("(1) DESIGNING COMMON METHYLATION CALLING WINDOWS FOR ADJACENT CLUSTERS")
  MethylationCallingWindows = Create_MethylationCallingWindows(TFBS_cluster_coordinates = TFBSs,
                                                               max_intercluster_distance = max_interTF_distance,
                                                               max_window_width = max_window_width,
                                                               min_cluster_width = min_cluster_width)
  message(paste0(length(MethylationCallingWindows), " METHYLATION CALLING WINDOWS DESIGNED"))
  
  message("(2) CALLING METHYLATION AND SORTING")
  parallel::mclapply(seq_along(MethylationCallingWindows), function(i){
    
    print(i)
    
    CurrentWindow = MethylationCallingWindows[i]
    ExperimentType = suppressMessages(SingleMoleculeFootprinting::DetectExperimentType(Samples = sample))
    
    CallContextMethylation(sampleSheet = sampleSheet,
                           sample = sample,
                           genome = genome,
                           RegionOfInterest = CurrentWindow,
                           coverage = coverage,
                           ConvRate.thr = ConvRate.thr,
                           returnSM = TRUE) -> Methylation
    
    if(length(Methylation[[1]]) == 0){
      return()
    }
    
    if (!is.null(CytosinesToMask)){
      
      message("Masking Cytosines")
      source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
      MaskSNPs2(Methylation = Methylation, 
                CytosinesToMaks = CytosinesToMask, 
                MaskSMmat = TRUE, 
                Experiment = ExperimentType) -> Methylation
      
    }
    
    Overlaps = findOverlaps(TFBSs, CurrentWindow)
    TFBSs_to_sort = TFBSs[queryHits(Overlaps)]
    
    lapply(seq_along(TFBSs_to_sort), function(j){
      if(ExperimentType == "NO"){
        Methylation[[2]] = lapply(Methylation[[2]], function(x){x$DGCHN})
      }
      
      SortReadsBySingleTF(MethSM = Methylation[[2]], 
                          TFBS = TFBSs_to_sort[j], 
                          bins = bins, 
                          coverage = sorting_coverage)
    }) -> SortedReads_window
    names(SortedReads_window) = names(TFBSs_to_sort)
    SortedReads_window
    
  }, mc.cores = cores) -> SortedReads
  SortedReads = unlist(SortedReads, recursive = FALSE)
  
  message("(3) CALCULATE STATE FREQUENCIES")
  Reduce(rbind,
         parallel::mclapply(seq_along(SortedReads), function(i){
           
           StateQuantification_tbl = StateQuantification(SortedReads = SortedReads[[i]], states = NULL)
           StateQuantification_tbl$TFBS_cluster = names(SortedReads[i])
           StateQuantification_tbl
           
         }, mc.cores = cores)) -> StateFrequency_tbl
  
  return(list(TFBSs, SortedReads, StateFrequency_tbl))
  
}



#' Convenience wrapper to sort single molecule according to TFBS clusters at multiple sites in the genome
#' 
#' The function starts from a list of single TFBSs, arranges them into clusters, calls methylation at the interested sites and outputs sorted reads
#' 
#' @param sampleSheet QuasR pointer file
#' @param sample samples to use, from the SampleName field of the sampleSheet
#' @param genome BSgenome
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to 0.8. To skip this filtering step, set to NULL. For more information, check out the details section.
# #' @param clObj cluster object for parallel processing of multiple samples/RegionsOfInterest. For now only used by qMeth call for bulk methylation. Should be the output of a parallel::makeCluster() call
#' @param CytosinesToMask THIS IS DUCK-TAPE
#' @param TFBSs GRanges object of transcription factor binding sites coordinates
#' @param max_intersite_distance maximum allowed distance in base pairs between two TFBS centers for them to be considered part of the same cluster. Defaults to 75.
#' @param min_intersite_distance minimum allowed distance in base pairs between two TFBS centers for them not to be discarded as overlapping. 
#'                               This parameter should be set according to the width of the bins used for later sorting. Defaults to 15.
#' @param max_cluster_size maximum number of TFBSs to be contained in any given cluster. Defaults to 10
#' @param max_intercluster_distance maximum distance between two consecutive TFBS clusters for them to be grouped in the same window
#' @param max_window_width upper limit to window width. This value should be adjusted according to the user's system as it determines the amount of memory used in the later context methylation call
#' @param min_cluster_width lower limit to window width. Corresponds to the scenario when a window contains a single TFBS cluster.
#' @param bins list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-7,7), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the most upstream TFBS.
#'             bins[[2]] represents all the TFBS bins, with coordinates relative to the center of each TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the most downstream TFBS.
#' @param sorting_coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed. Defaults to 30.
#' @param cores number of cores to use for parallel processing of multiple Methylation Calling Windows (i.e. groupings of adjecent TFBS clusters)
#' 
#' @importFrom parallel mclapply
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' 
#' @return list where [[1]] is the TFBS_Clusters object describing coordinates and composition of the TFBS clusters used to sort single molecules
#'                    [[2]] is a list of SortedReads nested per TFBS_cluster and sample
#'                    [[3]] is a tibble reporting the count (and frequency) of reads per state, sample and TFBS cluster
#' 
#' @export
#' 
SortReadsByTFCluster_MultiSiteWrapper = function(sampleSheet, sample, genome, coverage = 20, ConvRate.thr = 0.8, # clObj=NULL, ---> parameters passed to CallContextMethylation
                                                 CytosinesToMask = NULL,
                                                 TFBSs, max_intersite_distance = 75, min_intersite_distance = 15, max_cluster_size = 10,  # ---> parameters passed to Arrange_TFBSs_clusters
                                                 max_intercluster_distance = 100000, max_window_width = 5000000, min_cluster_width = 600, # ---> parameters passed to Create_MethylationCallingWindows
                                                 sorting_coverage = 30, bins = list(c(-35,-25), c(-7,7), c(25,35)), # ---> parameters passed to SortReadsByTFCluster
                                                 cores = 1
                                                 ){
  
  message("(1) CONSTRUCTING TFBS CLUSTERS")
  TFBS_Clusters = Arrange_TFBSs_clusters(TFBSs, 
                                         max_intersite_distance = max_intersite_distance,
                                         min_intersite_distance = min_intersite_distance, 
                                         max_cluster_size = max_cluster_size)
  message(paste0(length(TFBS_Clusters$ClusterCoordinates), " CLUSTERS FOUND"))
  
  message("(2) DESIGNING COMMON METHYLATION CALLING WINDOWS FOR ADJACENT CLUSTERS")
  MethylationCallingWindows = Create_MethylationCallingWindows(TFBS_cluster_coordinates = TFBS_Clusters$ClusterCoordinates,
                                                               max_intercluster_distance = max_intercluster_distance,
                                                               max_window_width = max_window_width,
                                                               min_cluster_width = min_cluster_width)
  message(paste0(length(MethylationCallingWindows), " METHYLATION CALLING WINDOWS DESIGNED"))
  
  message("(3) CALLING METHYLATION AND SORTING")
  parallel::mclapply(seq_along(MethylationCallingWindows), function(i){
    
    print(i)
    
    CurrentWindow = MethylationCallingWindows[i]
    ExperimentType = suppressMessages(SingleMoleculeFootprinting::DetectExperimentType(Samples = sample))
    
    CallContextMethylation(sampleSheet = sampleSheet,
                           sample = sample,
                           genome = genome,
                           RegionOfInterest = CurrentWindow,
                           coverage = coverage,
                           ConvRate.thr = ConvRate.thr,
                           returnSM = TRUE) -> Methylation
    
    if(length(Methylation[[1]]) == 0){
      return()
    }
    
    if (!is.null(CytosinesToMask)){
      
      message("Masking Cytosines")
      source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
      MaskSNPs2(Methylation = Methylation, 
               CytosinesToMaks = CytosinesToMask, 
               MaskSMmat = TRUE, 
               Experiment = ExperimentType) -> Methylation
      
    }
    
    Overlaps = findOverlaps(TFBS_Clusters$ClusterCoordinates, CurrentWindow)
    Clusters_to_sort = TFBS_Clusters$ClusterComposition[queryHits(Overlaps)]
    
    lapply(seq_along(Clusters_to_sort), function(j){
      if(ExperimentType == "NO"){
        Methylation[[2]] = lapply(Methylation[[2]], function(x){x$DGCHN})
      }
      SortReadsByTFCluster(MethSM = Methylation[[2]],
                           TFBS_cluster = Clusters_to_sort[[j]],
                           bins = bins, 
                           coverage = sorting_coverage)
    }) -> SortedReads_window
    names(SortedReads_window) = names(Clusters_to_sort)
    SortedReads_window
    
  }, mc.cores = cores) -> SortedReads
  SortedReads = unlist(SortedReads, recursive = FALSE)
  
  message("(4) CALCULATE STATE FREQUENCIES")
  Reduce(rbind,
         parallel::mclapply(seq_along(SortedReads), function(i){
           
           StateQuantification_tbl = StateQuantification(SortedReads = SortedReads[[i]], states = NULL)
           StateQuantification_tbl$TFBS_cluster = names(SortedReads[i])
           StateQuantification_tbl
           
         }, mc.cores = cores)) -> StateFrequency_tbl
  
  return(list(TFBS_Clusters, SortedReads, StateFrequency_tbl))
  
}

