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

  plot(NA,xlim=c(start(range),end(range)),ylim=c(-0.2,1),xlab='',ylab='SMF (%)',main=range)
  points(start(MethGR), 1-elementMetadata(MethGR)[[1]], type='l')
  points(start(MethGR), 1-elementMetadata(MethGR)[[1]], pch=20)
  abline(h=0)
  rect(start(TFBSs),-0.2,end(TFBSs),-0.15)
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
  plot(NA,pch='_',col=colors[as.factor(vR1[[3]])],xlim=c(start(range),end(range)),ylim=c(-1,length(unique(vR1[[2]]))))
  points(vR1[[1]],vR1[[2]],pch='_',cex=1.2,col=colors[as.factor(vR1[[3]])],xlim=c(start(range),end(range)))

}

#' Wrapper for PlotSingleMoleculeStack
#'
#' adds the covenience of arranging reads before plotting
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
    MethSM = MethSM[rev(unlist(SortedReads)),]
    PlotSingleMoleculeStack(MethSM, range)
  } else if (SortedReads == "HC"){
    message("Perfoming hierarchical clustering on single molecules before plotting")
    MethSM = HierarchicalClustering(MethSM)
    PlotSingleMoleculeStack(MethSM, range)
  }

}

#' Plot SMF data at single site
#'
#' @param ContextMethylation Context methylation object as returned by CallContextMethylation function
#' @param sample one sample as reported in the SampleName files of the QuasR sampleSheet
#' @param range GRange interval to plot
#' @param SortedReads Defaults to NULL, in which case will plot unsorted reads. Sorted reads object as returned by SortReads function or "HC" to perform hierarchical clustering
#' @param TFBSs GRange or GRangesList of transcription factor binding sites to add to the plot. If SortedReads are passed, the format of TFBSs (GRanges vs GRangesList) will be used to determie if single molecules were sorted based on one or multiple TFs
#' @param OutPath Full path to pdf file to save plot to. Defaults to NULL, in which case will only display
#'
#' @export
#'
PlotSingleTF = function(ContextMethylation, sample, range, SortedReads=NULL, TFBSs, OutPath=NULL){

  SortByCluster = tryCatch(length(TFBSs[[1]])>0, error = function(e){F})

  extende_range = resize(range,600,fix='center')
  # Subset TFBSs
  if (SortByCluster){
    subset_TFBSs = subsetByOverlaps(TFBSs, extende_range, ignore.strand=T)[[1]]
  } else {
    subset_TFBSs = subsetByOverlaps(TFBSs, extende_range, ignore.strand=T)
  }

  # Subset Cytosines
  MethGR = subsetByOverlaps(ContextMethylation[[1]], extende_range)
  MethSM = ContextMethylation[[2]][, as.character(start(ContextMethylation[[1]]))]

  ## PLOT ##
  # start graphical device
  if (!is.null(OutPath)){
    if(!is.null(dev.list())){dev.off()}
    pdf(OutPath, width = 8, height = 8)
  }

  # Average
  PlotAvgSMF(MethGR, extende_range, subset_TFBSs)

  # Single Molecule
  PlotSM(MethSM = MethSM, range = extende_range, SortedReads = SortedReads)

  if (is.list(SortedReads)){

    #  I GOT HERE <------------ NEED TO DO THE SORTING FIRST
    sRs = SortReads_internal(SortedReads, SM_mat = MethSM, isClusters = SortByCluster)
    readIDs = unlist(sRs, recursive = T, use.names = F)
    M=SM_mat[readIDs,]

    # Classification plot
    #add the classification plot

    if(SortByCluster){
      states = TFpairStates()

      TF1c=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(9)[c(9,9,4,9)]
      TF2c=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(9)[c(9,4,9,9)]
      names(TF1c)=c("00", "01", "10", "11")
      names(TF2c)=c("00", "01", "10", "11")

      boundaries=cumsum(lengths(sRs))
      TF1colv=lapply(seq(length(boundaries)),function(i){
        bd=c(0,boundaries)
        offset = as.integer(ifelse(bd[i+1] - bd[i] > 0, -1, 0)) # fixes errors for empty states
        rep(TF1c[substr(names(sRs),1,2)][i],length(seq(bd[i],bd[i+1]-offset,1)))
      })

      TF2colv=lapply(seq(length(boundaries)),function(i){
        bd=c(0,boundaries)
        offset = as.integer(ifelse(bd[i+1] - bd[i] > 0, -1, 0)) # fixes errors for empty states
        rep(TF2c[substr(names(sRs),3,4)][i],length(seq(bd[i],bd[i+1]-offset,1)))
      })

      plot(NA,xlim=c(0,3),ylim=c(0,sum(lengths(sRs))))
      points(rep(1, sum(lengths(sRs))-1),seq(sum(lengths(sRs))-1),col=(unlist(TF1colv)),pch='_',cex=2)
      points(rep(1.5, sum(lengths(sRs))-1),seq(sum(lengths(sRs))-1),col=(unlist(TF2colv)),pch='_',cex=2)

    } else {
      states = OneTFstates()

      counts=unlist(lapply(sRs,length))
      grouped.counts=unlist(lapply(seq_along(states),function(i){sum(counts[states[[i]]])}))
      names(grouped.counts)=names(states)

      TF1c=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(9)[c(2,3,4,9 )]
      names(TF1c)=names(states)
      boundaries=cumsum(grouped.counts)
      boundaries=boundaries
      TF1colv=lapply(seq(length(boundaries)),function(j){
        bd=grouped.counts
        rep(TF1c[names(grouped.counts)][j],bd[j])

      })

      plot(NA,xlim=c(0,3),ylim=c(0,sum(lengths(sRs))))
      points(rep(1, sum(grouped.counts)-1),seq(sum(grouped.counts)-1),col=unlist(TF1colv),pch='_',cex=2)
      text(rep(2,length(grouped.counts)),boundaries,round(grouped.counts/sum(grouped.counts)*100))
    }

  }

  if (!is.null(OutPath)){
    dev.off()
  }

}

