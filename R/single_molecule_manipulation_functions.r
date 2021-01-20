##################################################
#sort reads around TF binding sites
#to be used for factors binding in isolation
#requires presence of GC or CG in each of the three collection bins
#creates 3 bins
#output: list of sorted read IDs for each of the covered instance in the target_range (TFBS)
##################################################
sortTFreads_vect_m2<-function(methAln,selCs,st,inMs=c(NA,NA),upMs=c(NA,NA),doMs=c(NA,NA)){

	Crange=IRanges(methAln[[2]][selCs],methAln[[2]][selCs])
	elementMetadata(Crange)$meth = methAln[[4]][selCs]
	elementMetadata(Crange)$readID = methAln[[1]][selCs]

	#create collecting intervals
	binMethylationList = lapply(list(upMs, inMs, doMs), BinMethylation, Crange = Crange, st = st)

	# Subset the reads that cover all bins
	ReadsSubset = binMethylationList[[2]][[2]][binMethylationList[[2]][[2]] %in% Reduce(intersect, lapply(binMethylationList, function(x){x[[2]]}))]

	# binary met vectors
	binMethylationList_subset = lapply(binMethylationList, function(x){as.character(x[[1]][ReadsSubset])})
	pattern = Reduce(paste0, binMethylationList_subset)

	stID = unlist(lapply(strsplit(ReadsSubset, "_[A-Z]"), function(x){x[1]})) # If read ID doesn't start with capital letter this bugs
	ReadsID = unlist(lapply(strsplit(ReadsSubset, "_"), function(x){x[3]}))

	#SORTED READ LISTS
	if(length(ReadsSubset)>0){
	  sortedReads = split(ReadsID, paste(stID, pattern, sep='_'))
	  sortedReadslist = split(sortedReads, unlist(lapply(strsplit(names(sortedReads), "_"), function(x){paste(x[[1]], x[[2]], sep = "_")})))
	}else{
	  sortedReadslist=list()
	  }

	return(sortedReadslist)

  }

##################################################
#sort reads around clusters of TF binding sites
#to be used for factors binding in clusters
#creates n+2 bins (n=number of TFs in the cluster)
#output: list of sorted read IDs for each of the covered instance in the target_range (TFBS)
##################################################
sortTF_clustersreads_vect_m2<-function(methAln,selCs,st,inMs=c(NA,NA),upMs=c(NA,NA),doMs=c(NA,NA)){

  Crange=IRanges(methAln[[2]][selCs],methAln[[2]][selCs])

  # Make bins
  # in
  GRangesList(lapply(st, function(x){
    gr = resize(x, 1, fix='center')
    gr_resized = GRanges(seqnames(gr), IRanges(start = start(gr)+inMs[1], end = start(gr)+inMs[2]))
    names(gr_resized) = names(x)
    gr_resized
  })) -> inMPs
  # up
  upMPs=unlist(GRangesList(lapply(seq_along(st),function(i){
    TFcluster=st[[i]]
    minV=order(start(TFcluster),decreasing=F)[1]
    upMP=GRanges(seqnames(TFcluster[minV]),IRanges(start(TFcluster[minV])+upMs[1],start(TFcluster[minV])+upMs[2]))
    upMP
  })))
  names(upMPs)=names(st)
  # down
  doMPs=unlist(GRangesList(lapply(seq_along(st),function(i){
    TFcluster=st[[i]]
    maxV=order(start(TFcluster),decreasing=T)[1]
    doMP=GRanges(seqnames(TFcluster[maxV]),IRanges(start(TFcluster[maxV])+doMs[1],start(TFcluster[maxV])+doMs[2]))
    doMP
  })))
  names(doMPs)=names(st)

  #generate dynamic bins depending on the conposition of the TF cluster
  nb_TF=lengths(st)


  #group clusters by number of TFs
  #create dynamic bins
  #d_bins are grouped by binning type (number of TFs in the cluster)
  ClusterSizes = sort(unique(nb_TF))
  d_bins = lapply(seq_along(ClusterSizes),function(i){

    nbTF = ClusterSizes[i]
    TFrange = inMPs[lengths(inMPs) == nbTF]

    se = lapply(seq_along(TFrange), function(cl){
      lapply(seq(nbTF), function(pos){TFrange[[cl]][pos]})
    })

    inMPsd=lapply(seq(nbTF),function(pos){
    do.call(c,lapply(se, `[[`, pos))
    })

    cluster_id=names(st)[i]
    bins=GRangesList(c(list(upMPs[names(TFrange)]),inMPsd,list(doMPs[names(TFrange)])	))
    bins
  })
  names(d_bins) = paste0('clusterOf_', ClusterSizes, '_TF')


  #iterate ober bin types
  sRls=lapply(seq_along(d_bins),function(bin_t){

    #print(bin_t)
    bins=d_bins[[bin_t]]
    #find overlaping Cs
    ovs=lapply(seq_along(bins),function(i){
      as.matrix(findOverlaps(ranges(bins[[i]]),Crange)  )
    })
    #compute the methylation vectors
    Mvs=lapply(seq_along(bins),function(i){
      methAln[[4]][selCs][ovs[[i]][,2]]
    })
    #compute the readIDs
    IDvs=lapply(seq_along(bins),function(i){
      paste(names(bins[[1]])[ovs[[i]][,1]],methAln[[1]][selCs][ovs[[i]][,2]],sep='_')
    })
    #group Cs per region/per_read
    g.Ms=lapply(seq_along(bins),function(i){
      round(tapply(Mvs[[i]],IDvs[[i]],mean))
    })
    #find unique IDs
    u.IDs=lapply(seq_along(bins),function(i){
      sort(unique(IDvs[[i]]))
    })

    intMat=do.call(cbind,lapply(seq_along(bins),function(i){
      u.IDs[[1]]%in% u.IDs[[i]]
    }))

    #intersect IDS
    ids=seq_along( u.IDs[[1]])[rowSums(intMat)==ncol(intMat)] #makes sure the three bins are covered

    #  	s.outID=seq_along( u.outID)[u.outID %in% intersect( u.inID ,u.outID)]
    #binary met vectors
    sID=u.IDs[[1]][ids]


    s.M=lapply(seq_along(bins),function(i){
      g.Ms[[i]][sID]
    })

    patternMat=do.call(cbind,s.M)
    #print(patternMat)
    pattern=apply(patternMat,1,function(x){paste(as.character(x),sep='',collapse='')})

    st.id=paste(string.split(sID,'_',1),string.split(sID,'_',2),string.split(sID,'_',3),sep='_')
    read.id=string.split(sID,'_',4)
    #TABULAR COUNTS
    # 	p.counts=table(st.id,pattern)

    #SORTED READ LISTS
    if(length(ids)>0){
      sR=split(read.id,paste(st.id,pattern,sep='_'))
      sRl=split(sR,paste(string.split(names(sR),'_',1),string.split(names(sR),'_',2),string.split(names(sR),'_',3),sep='_'))

    }else{
      sRl = list()
    }

    sRl
  })
  #merge the results

  sRlsm = do.call(c, sRls)
}

#' Sort Single Molecules
#'
#' this function sorts molecules according to their state e.g. TF-bound, accessible, nucleosomee-occupied
#'
#' @export
#'
#' @param regDF Genomic chunk where single molecule methylation will be extracted (keep it ~1e8bp for optimal memory usage) provided as a data.frame
#' @param sampleSheet QuasR input file containing sample information (see \code{\link[QuasR]{qAlign}} documentation for more details). N.b. SampleNames have to be tagged with the type of treatment _NO_: only GCs will be used, _SS_: only CGs will be considered, _DE_: CG and GCs will be considered
#' @param target_range GRanges of TF binding motifs to be analyzed. N.b. If target range is a Grange, sorting will be carried in a single TF fashion. If is a GrangeList, reads will be sorted by TF clusters
#' @param inMs borders of the bin covering the TF binding site
#' @param upM borders of the bin upstream of the binding site (relative to the center of the motif)
#' @param doMs borders of the bin downstream of the binding site  (relative to the center of the motif)
#' @param conv.rate Integer from 0 to 100 representing the minimal accepted erroneous bisulfite conversion rate before discarding the read
#' @param remove.duplicates T/F. Remove identical reads, assumed to be PCR duplicates
SortSingleMolecules = function(regDF, sampleSheet, target_range, inMs=c(NA,NA), upMs=c(NA,NA), doMs=c(NA,NA), conv.rate=0, remove.duplicates=F){

  library("BSgenome.Mmusculus.UCSC.mm10")

  #load path to the alignments
  projAll <- QuasR::qAlign(sampleSheet,"BSgenome.Mmusculus.UCSC.mm10",paired="fr",bisulfite="dir")
  projAll@aligner = "Rbowtie"
  proj <- projAll[projAll@alignments$SampleName==regDF[1,4]]

  #define regions to extract methylation information
  reg <- GRanges(seqnames = Rle(regDF[1,1]),ranges = IRanges(regDF[1,2], end = regDF[1,3]),seqlengths=seqlengths(Mmusculus))
  regExp <- reg # expand the region at the end to catch all the pairs within the window
  end(regExp) <- pmin(end(regExp)+2000,seqlengths(regExp)[as.character(seqnames(regExp))])

  #call methylation vectors
  methAln <- QuasR::qMeth(proj,regExp,mode="allC",reportLevel="alignment")[[1]]

  if(length(methAln[[1]])>0){ #control that at leas one read is found in the region

    #filter PCR duplicates and low conversion reads
    if(!is.na(conv.rate)| remove.duplicates ){

      print(paste0("Filtering out reads with conversion rate > ", as.character(conv.rate)))
      if (remove.duplicates){print("Removing duplicates")}
      methAlnF = Filter_ConversionRate_and_Duplicates(methAln = methAln, reg = reg, conv.rate = conv.rate, remove.duplicates = remove.duplicates)
      names(methAlnF)=names(methAln)
      methAln=methAlnF
      print(paste0(as.character((1-(length(unique(methAlnF[[1]]))/length(unique(methAln[[1]]))))*100), "% reads filtered"))

    }

    selCs = CollapseStrands_SelectCytosines(methAln = methAln, regDF = regDF) # in the belly checks for context

    # Single TFs or clusters?
    SortByCluster = tryCatch(length(target_range[[1]])>0, error = function(e){F})
    if (SortByCluster){

      print("Sorting reads by TF clusters")

      #subset the target_range for ranges within the collection bin
      nb_ov = countOverlaps(unlist(target_range),GRanges(regDF),type='within')
      fact = as.factor(unlist(lapply(strsplit(names(unlist(target_range)), "_"), function(x){x[3]})))
      cluster.i=tapply(nb_ov,fact,sum)
      ids = paste0('TFBS_cluster_', unlist(lapply(strsplit(names(which(cluster.i>0)), "\\."), function(x){x[1]})))
      TFs=target_range[ids]

      if (any(is.na(c(inMs,upMs,doMs)))){
        print("Bins coordinates not passed, using default values for clusters: [-35,-25];[-7,7];[25,35]")
        inMs=c(-7,7)
        upMs=c(-35,-25)
        doMs=c(25,35)
      }

      sRs=sortTF_clustersreads_vect_m2(methAln,selCs,TFs,inMs,upMs,doMs)
      sRs

    } else {

      print("Sorting reads by single TFs")

      TFs=subsetByOverlaps(target_range,GRanges(regDF),type='within')

      if (any(is.na(c(inMs,upMs,doMs)))){
        print("Bins coordinates not passed, using default values for single TFs: [-35,-25];[-15,15];[25,35]")
        inMs=c(-15,15)
        upMs=c(-35,-25)
        doMs=c(25,35)
      }

      sRs=sortTFreads_vect_m2(methAln,selCs,TFs,inMs,upMs,doMs)
      sRs

    }

  }else{

    list()

  }

}
