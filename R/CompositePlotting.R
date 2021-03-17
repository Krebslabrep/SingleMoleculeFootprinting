#' Collect bulk SMF data for later composite plotting
#'
#' @param sampleSheet QuasR pointer file
#' @param sample for now this works for sure on one sample at the time only
#' @param genome BSgenome
#' @param TFBSs GRanges object of TF binding sites to collect info for
#' @param minCytosines minimum number of cytosines covered per TFBS
#' @param cores number of cores
#'
#' @import GenomicRanges
#' @import dplyr
#'
#' @return data.frame of bulk SMF info ready for plotting
#'
#' @export
#'
CollectCompositeData = function(sampleSheet, sample, genome, TFBSs, minCytosines=0, cores=20){

  # TO-DO: option to fix TFBS orientation | option to show C freq histogram | conditional plotting of geom_hex

  mclapply(seq_along(TFBSs), function(n){

    CallContextMethylation(sampleSheet = sampleSheet, sample = sample, genome = genome, range = TFBSs[n], returnSM = FALSE)

  }, mc.cores = cores) -> MethCalls
  names(MethCalls) = paste0("TFBS_", seq_along(TFBSs))

  GRlengths = lengths(MethCalls)
  FilteredMethCalls = MethCalls[GRlengths > minCytosines]
  message(paste0(round((length(FilteredMethCalls)/length(MethCalls))*100)), "% of sites are found covered at least at ", minCytosines, " cytosines")

  do.call(rbind, mclapply(seq_along(FilteredMethCalls), function(n){

    TFBS_index = names(FilteredMethCalls)[n]
    TFBS = TFBSs[as.integer(strsplit(TFBS_index, split = "_")[[1]][2])]
    TFBS_center = start(GenomicRanges::resize(TFBS, width = 1, fix = "center"))

    FilteredMethCalls[[n]] %>%
      as_tibble() %>%
      # select(-seqnames, -end, -width, -strand) %>%
      select(-end, -width, -strand) %>%
      mutate(relStart = start - TFBS_center, TFBS = TFBS_index)

  }, mc.cores = cores)) -> CompositeDF
  colnames(CompositeDF)[3] = "SMF"

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
