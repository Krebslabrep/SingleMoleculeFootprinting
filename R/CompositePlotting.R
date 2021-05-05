#' Collect bulk SMF data for later composite plotting
#'
#' @param sampleSheet QuasR pointer file
#' @param sample for now this works for sure on one sample at the time only
#' @param genome BSgenome
#' @param TFBSs GRanges object of TF binding sites to collect info for
#' @param window window size to collect methylation information for
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to 0.8. For more information, check out the details section
#' @param minCytosines minimum number of cytosines covered per TFBS
#' @param clObj cluster object for parallel processing. Should be the output of a parallel::makeCluster() call
#'
#' @import GenomicRanges
#' @import dplyr
#' @importFrom tidyr separate
#'
#' @return data.frame of bulk SMF info ready for plotting
#'
#' @export
#'
CollectCompositeData = function(sampleSheet, sample, genome, TFBSs, window, coverage=20, ConvRate.thr = 0.8, minCytosines=0, clObj=NULL){

  # TO-DO: option to fix TFBS orientation | option to show C freq histogram | conditional plotting of geom_hex
  
  TFBSs_extended = resize(TFBSs, width = window, 'center')

  MethCalls = CallContextMethylation(sampleSheet = sampleSheet, 
                                     sample = sample, 
                                     genome = genome, 
                                     RegionOfInterest = TFBSs_extended, 
                                     coverage = coverage, 
                                     ConvRate.thr = ConvRate.thr,
                                     returnSM = FALSE, 
                                     clObj = clObj)
  Overlaps = findOverlaps(MethCalls, TFBSs)
  MethCalls$TFBS_index = NA
  MethCalls[queryHits(Overlaps)]$TFBS_index = subjectHits(Overlaps)

  message(paste0(round(sum(table(MethCalls$TFBS_index) > 3)/length(TFBSs)*100), 
                 "% of sites are found covered at least at ", minCytosines, " cytosines"))
  
  MethCalls %>%
    as_tibble() %>%
    select(-end, -width) %>%
    mutate(TFBS_center = start(resize(TFBSs[subjectHits(Overlaps)], width = 1, fix = 'center'))) %>%
    mutate(RelStart = start - TFBS_center) %>% 
    gather(Measure, Value, -seqnames, -start, -strand, -GenomicContext, -TFBS_index, -TFBS_center, -RelStart) %>%
    mutate(Sample = gsub("_MethRate$|_Coverage$", "", Measure), Measure = gsub("^.*_", "", Measure)) %>%
    spread(Measure, Value) -> CompositeDF
  
  
  # MethCalls
  # 
  # do.call(rbind, mclapply(seq_along(FilteredMethCalls), function(n){
  # 
  #   TFBS_index = names(FilteredMethCalls)[n]
  #   TFBS = TFBSs[as.integer(strsplit(TFBS_index, split = "_")[[1]][2])]
  #   TFBS_center = start(GenomicRanges::resize(TFBS, width = 1, fix = "center")))
  # 
  #   FilteredMethCalls[[n]] %>%
  #     as_tibble() %>%
  #     # select(-seqnames, -end, -width, -strand) %>%
  #     select(-end, -width, -strand) %>%
  #     mutate(relStart = start - TFBS_center, TFBS = TFBS_index)
  # 
  # }, mc.cores = cores)) -> CompositeDF
  # colnames(CompositeDF)[3] = "SMF"

  return(CompositeDF)

}
#
# n = length(top500scores)
# TF = as.character(unique(top500scores$name))
# CompositeDF = CollectCompositeData(sampleSheet = Qinput, sample = MySample[1], genome = BSgenome.Mmusculus.UCSC.mm10, TFBSs = top500scores, minCytosines = 5, cores = 50)
#
# library(ggplot2)
# library(viridis)
# CompositeDF %>%
#   ggplot() +
#   geom_hex(aes(start, SMF)) +
#   geom_smooth(aes(start, SMF), se = TRUE, colour = "#00BBDB") +
#   scale_fill_viridis(option = "inferno") +
#   ylim(c(0,1)) +
#   xlab("Coord relative to TFBS center") +
#   ggtitle(paste0("Top ", n, " bound ", TF)) +
#   theme_classic() +
#   ggsave("/g/krebs/barzaghi/analyses/misc/TP53_check/TP53_composite.pdf", width = 7, height = 5)
