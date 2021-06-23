#' Collect bulk SMF data for later composite plotting
#'
#' @param sampleSheet QuasR pointer file
#' @param sample for now this works for sure on one sample at the time only
#' @param genome BSgenome
#' @param TFBSs GRanges object of TF binding sites to collect info for. We reccommend employing 50 to 200 TFBSs.
#' @param window window size to collect methylation information for
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to 0.8. For more information, check out the details section
#' @param minCytosines minimum number of cytosines covered per TFBS
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
CollectCompositeData = function(sampleSheet, sample, genome, TFBSs, window, coverage=20, ConvRate.thr = 0.8, minCytosines=0, cores=1){

  # TO-DO: option to fix TFBS orientation | option to show C freq histogram | conditional plotting of geom_hex
  
  TFBSs_extended = resize(TFBSs, width = window, 'center')
  
  Reduce(c,
  parallel::mclapply(seq_along(TFBSs_extended), function(i){

    CallContextMethylation(sampleSheet = sampleSheet, 
                           sample = sample, 
                           genome = genome, 
                           RegionOfInterest = TFBSs_extended[i], 
                           coverage = coverage, 
                           ConvRate.thr = ConvRate.thr,
                           returnSM = FALSE, 
                           clObj = NULL) -> Methylation
    if(length(Methylation)>0){
      Methylation$TFBS_index = i
      Methylation$TFBS_center = start(resize(TFBSs_extended[i], width = 1, fix = 'center'))
    } else {
      Methylation$TFBS_index = c()
      Methylation$TFBS_center = c()
    }
    Methylation
    
  }, mc.cores = cores)) -> MethCalls

  message(paste0(round(sum(table(MethCalls$TFBS_index) > minCytosines)/length(TFBSs)*100), 
                 "% of sites are found covered at least at ", minCytosines, " cytosines"))
  
  MethCalls %>%
    as_tibble() %>%
    mutate(RelStart = start - TFBS_center) %>% 
    select(-seqnames, -start, -end, -strand, -width, -TFBS_center, -GenomicContext) %>%
    gather(Measure, Value, -RelStart, -TFBS_index) %>%
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







