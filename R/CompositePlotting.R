#' Collect bulk SMF data for later composite plotting
#'
#' @param sampleSheet QuasR pointer file
#' @param sample for now this works for sure on one sample at the time only
#' @param genome BSgenome
#' @param TFBSs GRanges object of TF binding sites to collect info for. We reccommend employing 50 to 200 TFBSs.
#' @param window window size to collect methylation information for
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to NULL For more information, check out the details section
#' @param minCytosines minimum number of cytosines covered per TFBS
#' @param max_window_width Defautls to 5e+06. Check Create_MethylationCallingWindows documentation for details
#' @param max_intercluster_distance Defautls to 1e+06. Check Create_MethylationCallingWindows documentation for details
#' @param cores number of cores to use
#'
#' @import GenomicRanges
#' @import dplyr
#' @importFrom tidyr separate
#'
#' @return data.frame of bulk SMF info ready for plotting
#'
#' @export
#'
CollectCompositeData = function(sampleSheet, sample, genome, TFBSs, window, coverage=20, ConvRate.thr = NULL, minCytosines=0, max_window_width = 5e+06, max_intercluster_distance = 1e+06, cores=1){

  # TO-DO: option to fix TFBS orientation | option to show C freq histogram | conditional plotting of geom_hex
  
  TFBSs_extended = resize(TFBSs, width = window, 'center')
  
  MethylationCallingWindows = Create_MethylationCallingWindows(TFBS_cluster_coordinates = TFBSs_extended, max_window_width = max_window_width, max_intercluster_distance = max_intercluster_distance, genomic.seqlenghts = GenomeInfoDb::seqlengths(genome))
  
  Reduce(c,
  parallel::mclapply(seq_along(MethylationCallingWindows), function(i){
    
    CallContextMethylation(sampleSheet = sampleSheet, 
                           sample = sample, 
                           genome = genome, 
                           RegionOfInterest = MethylationCallingWindows[i], 
                           coverage = coverage, 
                           ConvRate.thr = ConvRate.thr,
                           returnSM = FALSE, 
                           clObj = NULL) -> Methylation
    Methylation$TFBS_center = NA
    
    findOverlaps(Methylation, TFBSs_extended) -> Overlaps
    if(length(Overlaps) == 0){
      return(Methylation)
    }
    
    # This also expands cytosines overlapping with more than 1 TFBS    
    Methylation.of.interest = Methylation[queryHits(Overlaps)]
    Methylation.of.interest$TFBS_center = start(IRanges::resize(TFBSs_extended[subjectHits(Overlaps)], width = 1, fix = 'center'))
    Methylation.of.interest$TFBS_index = subjectHits(Overlaps) #paste(seqnames(Methylation.of.interest), Methylation.of.interest$TFBS_center, sep = "_")
    
    return(Methylation.of.interest)
    
  }, mc.cores = cores, mc.preschedule = FALSE)) -> MethCalls
    
  message(paste0(round(sum(table(MethCalls$TFBS_index) > minCytosines)/length(TFBSs)*100), 
                 "% of sites are found covered at least at ", minCytosines, " cytosines"))
  
  MethCalls %>%
    as_tibble() %>%
    mutate(RelStart = start - TFBS_center) %>% 
    gather(Measure, Value, grep("_Coverage$|_MethRate$", colnames(elementMetadata(MethCalls)), value = TRUE)) %>%
    mutate(Sample = gsub("_MethRate$|_Coverage$", "", Measure), Measure = gsub("^.*_", "", Measure)) %>%
    spread(Measure, Value) %>%
    mutate(SMF = 1 - MethRate) %>%
    select(-MethRate) -> CompositeDF
  
  return(CompositeDF)

}


#' Plot composite SMF data
#'
#' Will use geom_point with <= 5000 points, geom_hex otherwise 
#' 
#' @param CompositeDF the output of the CollectCompositeData function
#' @param span the span parameter to pass to geom_smooth
#' 
#' @import ggplot2
#' @import viridis
#' 
#' @export
#' 
CompositePlot = function(CompositeDF, span=0.1, TF){
  
  if(nrow(CompositeDF) <= 5000){
    PlotHex = FALSE
  } else {
    PlotHex = TRUE
  }
  
  CompositeDF %>%
    ggplot() +
    {if(!PlotHex){geom_point(aes(RelStart, SMF), alpha=0.5)}} +
    {if(PlotHex){geom_hex(aes(RelStart, SMF))}} +
    geom_smooth(aes(RelStart, SMF), se = TRUE, colour = "#00BBDB", method = "loess", span = span) +
    scale_fill_viridis(option = "inferno") +
    ylim(c(0,1)) +
    xlab("Coord relative to TFBS center") +
    ggtitle(paste0("Top ", length(unique(CompositeDF$TFBS_index)), " covered ", TF, " sites")) +
    theme_classic() 
  
}







