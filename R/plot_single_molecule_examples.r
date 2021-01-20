# one sample, one TF/cluster
PlotSingleTF = function(AvgMeth, SingleMoleculeMeth, sample, range, TFBSs, SortedReads=NULL, OutPath=NULL){

  SortByCluster = tryCatch(length(TFBSs[[1]])>0, error = function(e){F})

  extende_range = resize(range,600,fix='center')
  # Subset TFBSs
  if (SortByCluster){
    subset_TFBSs = subsetByOverlaps(TFBSs, extende_range, ignore.strand=T)[[1]]
  } else {
    subset_TFBSs = subsetByOverlaps(TFBSs, extende_range, ignore.strand=T)
  }
  # Subset Cytosines
  AvgMeth_gr = subsetByOverlaps(AvgMeth, extende_range)
  SM_mat = SingleMoleculeMeth[, as.character(start(AvgMeth_gr))]

  ## PLOT ##
  # start graphical device
  if (!is.null(OutPath)){
    if(!is.null(dev.list())){dev.off()}
    pdf(OutPath, width = 8, height = 8)
  }

  # Average
  plot(NA,xlim=c(start(extende_range),end(extende_range)),ylim=c(-0.2,1),xlab='',ylab='SMF (%)',main=extende_range)
  points(colnames(SM_mat), 1-colMeans(SM_mat, na.rm='T'), type='l')
  points(colnames(SM_mat), 1-colMeans(SM_mat, na.rm='T'), pch=20)
  abline(h=0)
  rect(start(subset_TFBSs),-0.2,end(subset_TFBSs),-0.15)
  text(start(resize(subset_TFBSs,1,fix='center')),rep(-0.1,length(subset_TFBSs)),subset_TFBSs$name,cex=0.8)

  # Single Molecules stack
  sRs = SortReads_internal(SortedReads, SM_mat = SM_mat, isClusters = SortByCluster)
  readIDs = unlist(sRs, recursive = T, use.names = F)

  # Use Hierarchical Clustering if no previously sorted reads were passed
  if (is.null(SortedReads)){

    print("No sorted reads passed, proceeding with hierarcical clustering")
    if(length(readIDs)>500){ #subset to 500 molecules to avoid problem with Hc
      SM_mat_subset = SM_mat[sample(readIDs,500),]
    }else{
      SM_mat_subset = SM_mat[readIDs,]
    }
    hc=hclust(dist(SM_mat_subset))
    M=SM_mat_subset[hc$order,]

  } else {

    M=SM_mat[readIDs,]

  }

  vR1=VectorizeReads(extende_range,M)
  BR=c(col=c('black','grey'))
  colors=BR
  plot(NA,pch='_',col=colors[as.factor(vR1[[3]])],xlim=c(start(extende_range),end(extende_range)),ylim=c(-1,length(unique(vR1[[2]]))))
  points(vR1[[1]],vR1[[2]],pch='_',cex=1.2,col=colors[as.factor(vR1[[3]])],xlim=c(start(extende_range),end(extende_range)))

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


  if (!is.null(OutPath)){
    dev.off()
  }

}

