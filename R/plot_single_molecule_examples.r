#' Perform Hierarchical clustering on single reads
#'
#' @param MethSM Single molecule methylation matrix
#'
#' @importFrom stats hclust
#' @importFrom stats dist
#'
HierarchicalClustering = function(MethSM){

  if(nrow(MethSM) > 500){ #subset to 500 molecules to avoid problem with Hc
    MethSM_subset = MethSM[sample(dimnames(MethSM)[[1]],500),]
  }else{
    MethSM_subset = MethSM
  }
  hc=hclust(dist(MethSM_subset))
  MethSM_HC = MethSM_subset[hc$order,]

  return(MethSM_HC)

}

#' Plot average methylation
#'
#' @param MethGR Average methylation GRanges obj
#' @param range GRanges interval to plot
#' @param TFBSs GRanges object of transcription factor binding sites to include in the plot. Assumed to be already subset.
#'
#' @import GenomicRanges
#'
#' @export
#'
PlotAvgSMF = function(MethGR, range, TFBSs){

  plot(NA,xlim=c(start(range),end(range)),ylim=c(-0.2,1),xlab='',ylab='SMF',main=range)
  points(start(MethGR), 1-elementMetadata(MethGR)[[1]], type='l') #, lwd=4
  points(start(MethGR), 1-elementMetadata(MethGR)[[1]], pch=20) #, lwd=5
  abline(h=0)
  rect(start(TFBSs),-0.2,end(TFBSs),-0.15)#, lwd=2
  text(start(resize(TFBSs,1,fix='center')),rep(-0.1,length(TFBSs)),TFBSs$name,cex=0.8)

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
PlotSM = function(MethSM, range, SortedReads = NULL){

  if (is.null(SortedReads)){
    message("No sorting passed or specified, will plot unsorted reads")
    PlotSingleMoleculeStack(MethSM, range)
  } else if (is.list(SortedReads)){
    message("Sorting reads according to passed values before plotting")
    if (length(SortedReads) <= 8){ # Single TF
      message("Inferring sorting was performed by single TF")
      OrderedReads = SortedReads[as.character(unlist(OneTFstates()))]
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
#' @export
#'
SingleTFStateQuantificationPlot = function(states, OrderedReads){

  GroupedCounts = unlist(lapply(seq_along(states), function(i){
    length(unlist(OrderedReads[states[[i]]]))
  }))
  names(GroupedCounts) = names(states)

  Colors = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(9)[c(2,3,4,9 )]
  names(Colors) = names(states)
  boundaries = cumsum(GroupedCounts)
  ColorVector = lapply(seq_along(boundaries), function(j){
    rep(Colors[names(GroupedCounts)][j],GroupedCounts[j])
  })

  plot(NA,xlim=c(0,3),ylim=c(0,sum(lengths(OrderedReads))), ylab=paste0(nrow(MethSM), " reads"), xlab="")
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

