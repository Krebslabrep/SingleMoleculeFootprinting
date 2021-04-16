#' Perform Hierarchical clustering on single reads
#'
#' @param MethSM Single molecule methylation matrix
#'
#' @importFrom stats hclust
#' @importFrom stats dist
#'
HierarchicalClustering = function(MethSM){

  SubsetSize = 500
  if(nrow(MethSM) > SubsetSize){ #subset to 500 molecules to avoid problem with Hc
    MethSM_subset = MethSM[sample(dimnames(MethSM)[[1]],SubsetSize),]
  }else{
    SubsetSize = nrow(MethSM)
    MethSM_subset = MethSM
  }
  ReadsDist = dist(MethSM_subset)
  iteration = 0
  while(sum(is.na(ReadsDist)) > 0){ # sometimes dist between some pair of reads is NA, possibly because of no overlapping Cs
    iteration = iteration + 1
    if (iteration > SubsetSize) {
      SubsetSize = ceiling(SubsetSize*0.9)
      iteration = 0
    }
    MethSM_subset = MethSM[sample(dimnames(MethSM)[[1]],SubsetSize),]
    ReadsDist = dist(MethSM_subset)
  }
  hc=hclust(ReadsDist)
  MethSM_HC = MethSM_subset[hc$order,]

  return(MethSM_HC)

}

#' Plot average methylation
#'
#' @param MethGR Average methylation GRanges obj
#' @param RegionOfInterest GRanges interval to plot
#' @param TFBSs GRanges object of transcription factor binding sites to include in the plot. Assumed to be already subset. Also assumed that the tf names are under the column "TF"
#' @param SNPs GRanges object of SNPs to visualize. Assumed to be already subset. Assumed to have the reference and alternative sequences respectively under the columns "R" and "A"
#' @param SortingBins GRanges object of sorting bins (absolute) coordinate to visualize
#'
#' @import GenomicRanges
#' @import tidyverse
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#'
#' PlotAvgSMF(MethGR = Methylation[[1]], Region_of_interest = Region_of_interest, TFBSs = TFBSs)
#'
PlotAvgSMF = function(MethGR, RegionOfInterest, TFBSs=NULL, SNPs=NULL, SortingBins=NULL){

  # Prepare SMF data
  MethGR %>%
    as_tibble() %>%
    select(-grep("_Coverage$", colnames(.)), -end, -width, -strand) %>%
    gather(sample, MethRate, -seqnames, -start, -GenomicContext) -> PlottingDF
  OurFavouriteColors = c("Black", "Red", "Blue", "Green")
  ColorsToUse = OurFavouriteColors[seq_along(unique(PlottingDF$sample))]

  # Prepare TFBS
  if(!is.null(TFBSs)){
    TFBSs %>%
      as_tibble() %>%
      select(start, end, TF) -> TFBS_PlottingDF
  }

  # Prepare SNPs
  if(!is.null(SNPs)){
    SNPs %>%
      as_tibble() %>%
      select(start, R, A) %>%
      gather(Genotype, Sequence, -start) %>%
      mutate(y_coord = rep(c(-0.18,-0.21), each=length(SNPs))) %>%
      add_row(start=min(start(SNPs))-40, Genotype = c("R", "A"),
              Sequence =  c("Genotype:R", "Genotype:A"), y_coord = c(-0.18, -0.21)) -> SNPs_PlottingDF
  }

  # Prepare SortingBins
  if(!is.null(SortingBins)){
    SortingBins %>%
      as_tibble() %>%
      select(start, end) -> Bins_PlottingDF
  }

  PlottingDF %>%
    ggplot(aes(x=start, y=1-MethRate, color=sample)) +
    geom_line() +
    geom_point() +
    {if(!is.null(TFBSs)){geom_rect(TFBS_PlottingDF, mapping = aes(xmin=start, xmax=end, ymin=-0.15, ymax=-0.1), inherit.aes = FALSE)}} +
    {if(!is.null(TFBSs)){geom_text(TFBS_PlottingDF, mapping = aes(x=start+((end-start)/2), y=-0.08, label=TF), inherit.aes = FALSE)}} +
    {if(!is.null(SNPs)){geom_text(SNPs_PlottingDF, mapping = aes(x=start, y=y_coord, label=Sequence), size=3, inherit.aes = FALSE)}} +
    {if(!is.null(SortingBins)){geom_rect(Bins_PlottingDF, mapping = aes(xmin=start, xmax=end, ymin=-0.02, ymax=0), color="black", fill="white", inherit.aes = FALSE)}} +
    geom_hline(aes(yintercept=0)) +
    ylab("SMF") +
    xlab("") +
    ylim(c(-0.25,1)) +
    ggtitle(RegionOfInterest) +
    scale_color_manual(values = ColorsToUse) +
    theme_classic()

  # plot(NA,xlim=c(start(range),end(range)),ylim=c(-0.2,1),xlab='',ylab='SMF',main=range)
  # points(start(MethGR), 1-elementMetadata(MethGR)[[1]], type='l') #, lwd=4
  # points(start(MethGR), 1-elementMetadata(MethGR)[[1]], pch=20) #, lwd=5
  # abline(h=0)
  # rect(start(TFBSs),-0.2,end(TFBSs),-0.15)#, lwd=2
  # text(start(resize(TFBSs,1,fix='center')),rep(-0.1,length(TFBSs)),TFBSs$name,cex=0.8)

}

#' Plot single molecule stack
#'
#' @param MethSM Single molecule methylation matrix
#' @param range GRanges interval to plot
#'
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#'
#' PlotSingleMoleculeStack(MethSM = Methylation[[2]], range = Region_of_interest)
#'
PlotSingleMoleculeStack = function(MethSM, range){

  vR1=VectorizeReads(range,MethSM)
  BR=c(col=c('black','grey'))
  colors=BR
  plot(NA,pch='_',col=colors[as.factor(vR1[[3]])],xlim=c(start(range),end(range)),ylim=c(-1,length(unique(vR1[[2]]))), ylab=paste0(nrow(MethSM), " reads"), xlab="")
  points(vR1[[1]],vR1[[2]],pch='_',cex=1.2,col=colors[as.factor(vR1[[3]])],xlim=c(start(range),end(range)))

}

#' Wrapper for PlotSingleMoleculeStack function
#'
#' adds the convenience of arranging reads before plotting
#'
#' @param MethSM Single molecule methylation matrix
#' @param range GRanges interval to plot
#' @param SortedReads Defaults to NULL, in which case will plot unsorted reads. Sorted reads object as returned by SortReads function or "HC" to perform hierarchical clustering
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#'
#'  PlotSM(MethSM = Methylation[[2]], range = Region_of_interest)
#'
PlotSM = function(MethSM, range, SortedReads = NULL){

  if (is.null(SortedReads)){
    message("No sorting passed or specified, will plot unsorted reads")
    PlotSingleMoleculeStack(MethSM, range)
  } else if (is.list(SortedReads)){
    message("Sorting reads according to passed values before plotting")
    if (length(SortedReads) <= 8){ # Single TF
      message("Inferring sorting was performed by single TF")
      OrderedReads = SortedReads[as.character(unlist(OneTFstates()))]
      OrderedReads = rev(OrderedReads)
    } else if (length(SortedReads) > 8 & length(SortedReads) <= 16){ # TF pair
      message("Inferring sorting was performed by TF pair")
      OrderedReads = SortedReads[as.character(unlist(TFpairStates()))]
    }
    MethSM = MethSM[unlist(OrderedReads),]
    PlotSingleMoleculeStack(MethSM, range)
  } else if (SortedReads == "HC"){
    message("Perfoming hierarchical clustering on single molecules before plotting")
    MethSM = HierarchicalClustering(MethSM)
    PlotSingleMoleculeStack(MethSM, range)
  }

}

#' Single TF state quantification bar
#'
#' @param states as returned by OneTFstates function
#' @param OrderedReads Reads ordered by states
#'
#' @importFrom RColorBrewer brewer.pal
#'
SingleTFStateQuantificationPlot = function(states, OrderedReads){

  GroupedCounts = unlist(lapply(seq_along(states), function(i){
    length(unlist(OrderedReads[states[[i]]]))
  }))
  names(GroupedCounts) = names(states)

  Colors = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(9)[c(4,3,2,9)]
  names(Colors) = names(states)
  boundaries = cumsum(GroupedCounts)
  ColorVector = lapply(seq_along(boundaries), function(j){
    rep(Colors[names(GroupedCounts)][j],GroupedCounts[j])
  })

  GroupedCounts = rev(GroupedCounts)
  ColorVector = rev(ColorVector)
  boundaries = rev(boundaries)

  plot(NA,xlim=c(0,3),ylim=c(0,sum(GroupedCounts)), ylab=paste0(sum(GroupedCounts), " reads"), xlab="")
  points(rep(1, sum(GroupedCounts)-1),seq(sum(GroupedCounts)-1),col=unlist(ColorVector),pch='_',cex=2)
  text(rep(2,length(GroupedCounts)),boundaries,round(GroupedCounts/sum(GroupedCounts)*100))

}

#' TF pair state quantification bar
#'
#' @param states as returned by TFpairStates function
#' @param OrderedReads Reads ordered by states
#'
#' @importFrom RColorBrewer brewer.pal
#'
TFPairStateQuantificationPlot = function(states, OrderedReads){

  OrderedReads = OrderedReads[lengths(OrderedReads) > 0]

  GroupedCounts = unlist(lapply(seq_along(states), function(i){
    length(unlist(OrderedReads[states[[i]]]))
  }))
  names(GroupedCounts) = names(states)

  TF1c = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(9)[c(9,9,4,9)]
  TF2c = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(9)[c(9,4,9,9)]
  names(TF1c) = c("00", "01", "10", "11")
  names(TF2c) = c("00", "01", "10", "11")

  boundaries = cumsum(lengths(OrderedReads))
  TF1colv = lapply(seq_along(boundaries),function(i){
    bd = c(0,boundaries)
    rep(TF1c[substr(names(OrderedReads),1,2)][i],length(seq(bd[i],bd[i+1]-1,1)))
  })

  TF2colv = lapply(seq_along(boundaries),function(i){
    bd = c(0,boundaries)
    rep(TF2c[substr(names(OrderedReads),3,4)][i],length(seq(bd[i],bd[i+1]-1,1)))
  })

  TotReads = sum(lengths(OrderedReads))

  plot(NA, xlim = c(0,3), ylim = c(0,TotReads), ylab=paste0(TotReads, " reads"), xlab="")
  points(rep(1, TotReads-1), seq(TotReads-1), col=(unlist(TF1colv)), pch='_', cex=2)
  points(rep(1.5, TotReads-1), seq(TotReads-1), col=(unlist(TF2colv)), pch='_', cex=2)
  text(rep(2,length(GroupedCounts)),cumsum(GroupedCounts),round(GroupedCounts/sum(GroupedCounts)*100))

}

#' Plot states quantification bar
#'
#' @param SortedReads Sorted reads object as returned by SortReads function
#'
#' @return Bar plot quantifying states
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#' SortedReads = SortReadsByTFCluster(MethSM = Methylation[[2]], TFBSs = TFBSs)
#'
#' StateQuantificationPlot(SortedReads = SortedReads)
#'
StateQuantificationPlot = function(SortedReads){

  if (length(SortedReads) <= 8){ # Single TF

    message("Inferring sorting was performed by single TF")
    states = OneTFstates()
    OrderedReads = SortedReads[as.character(unlist(states))]
    SingleTFStateQuantificationPlot(states, OrderedReads)

  } else if (length(SortedReads) > 8 & length(SortedReads) <= 16){ # TF pair

    message("Inferring sorting was performed by TF pair")
    states = TFpairStates()
    OrderedReads = SortedReads[as.character(unlist(states))]
    TFPairStateQuantificationPlot(states, OrderedReads)

  }

}

#' Plot SMF data at single site
#'
#' @param ContextMethylation Context methylation object as returned by CallContextMethylation function
#' @param sample one sample as reported in the SampleName files of the QuasR sampleSheet
#' @param range GRange interval to plot
#' @param SortedReads Defaults to NULL, in which case will plot unsorted reads. Sorted reads object as returned by SortReads function or "HC" to perform hierarchical clustering
#' @param TFBSs GRange or GRangesList of transcription factor binding sites to add to the plot. If SortedReads are passed, the format of TFBSs (GRanges vs GRangesList) will be used to determie if single molecules were sorted based on one or multiple TFs
#' @param saveAs Full path to pdf file to save plot to. Defaults to NULL, in which case will only display
#'
#' @importFrom IRanges subsetByOverlaps resize
#' @importFrom GenomicRanges start
#' @import grDevices
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#' SortedReads = SortReadsByTFCluster(MethSM = Methylation[[2]], TFBSs = TFBSs)
#'
#' PlotSingleSiteSMF(ContextMethylation = Methylation,
#'                   sample = MySample,
#'                   range = Region_of_interest,
#'                   SortedReads = SortedReads,
#'                   TFBSs = TFBSs,
#'                   saveAs = NULL)
#'
PlotSingleSiteSMF = function(ContextMethylation, sample, range, SortedReads=NULL, TFBSs, saveAs=NULL){

  extende_range = resize(range, 600, fix='center')

  message("Subsetting data by range (extended)")
  subset_TFBSs = subsetByOverlaps(TFBSs, extende_range, ignore.strand=TRUE)
  MethGR = subsetByOverlaps(ContextMethylation[[1]], extende_range)
  MethSM = ContextMethylation[[2]][, as.character(start(ContextMethylation[[1]]))]

  ## PLOT ##
  # start graphical device
  if (!is.null(saveAs)){
    if(!is.null(dev.list())){dev.off()}
    pdf(saveAs, width = 8, height = 5)
  }

  # Average
  PlotAvgSMF(MethGR, extende_range, subset_TFBSs)

  # Single Molecule
  PlotSM(MethSM = MethSM, range = extende_range, SortedReads = SortedReads)

  # State quantification plot
  StateQuantificationPlot(SortedReads = SortedReads)

  if (!is.null(saveAs)){
    dev.off()
  }

}

