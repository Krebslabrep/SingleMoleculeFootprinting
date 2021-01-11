getCMethMatrix<-function(proj,range,samp){
	Cs=QuasR::qMeth(proj, query=range,mode="allC",reportLevel="alignment")
	# use data.table to get a 1,0 matrix of methylation profiles
	all.cids=unique(Cs[[samp]]$Cid) # get all possible C locations
	# make the data.table object
	dt=data.table(meth=Cs[[samp]]$meth ,aid=Cs[[samp]]$aid ,cid=Cs[[samp]]$Cid)
	# this function converts cids to columns
	myfun2<-function(x,all.cids){
	vec=rep(-1,length(all.cids))
	names(vec)=as.character(all.cids)
	b=as.list((vec))
	b[ as.character(x$cid)]=as.double(x$meth)
	return(b)
}
	dtm=dt[,myfun2(.SD,all.cids), by=aid]
	ronames=dtm$aid
	dtm[,aid:=NULL] # remove unwanted row
	CpGm=as.matrix(dtm)
	CpGm[CpGm == -1]=NA # put NAs
	rownames(CpGm)=ronames
	return(CpGm)
}

mergeGC_CGmat=function(GC_CG_mat){
			CGmat=GC_CG_mat$matCG
			GCmat=GC_CG_mat$matGC
			uReads=unique(c(rownames(CGmat),rownames(GCmat)))
			uCs=sort(unique(c(colnames(CGmat),colnames(GCmat))))
			matCGo=matrix(nrow=length(uReads),ncol=length(uCs))
			colnames(matCGo)=uCs
			rownames(matCGo)=uReads
			matCGo[rownames(CGmat),colnames(CGmat)]=CGmat#get CGs in
			matCGo[rownames(GCmat),colnames(GCmat)]=GCmat#get GCs in
	 		matCGo

		}

findContextCytosines <- function(meth_gr, genome, context){

  chr = as.character(unique(seqnames(meth_gr)))
  GenomicCytosines = matchPattern(context, genome[[chr]])
  GenomicCytosines_gr = GRanges(seqnames = chr, ranges = GenomicCytosines)
  meth_gr[meth_gr %over% GenomicCytosines_gr]

}

CollapeStrands <- function(meth_gr, context){

  # Fix Granges for overhangs: CG --> if first range is "-" add "+" with zeros on top, if it ends in "+" add "-" with zeros at the bottom
  if (length(meth_gr) > 0){ # check if the obj is empty

    TopStrandToFix = ifelse(context == "CG", "-", "+")
    if (as.character(strand(meth_gr[1]))==TopStrandToFix){
      print("Strand collapsing: Fixing top (left) overhang")
      overhang_fix = GRanges(seqnames = unique(as.character(seqnames(meth_gr))), ranges = IRanges(start(meth_gr[1])-1, width = 1), strand = ifelse(context == "CG", "+", "-"))
      fix_metadata = data.frame(matrix(0, ncol = length(colnames(elementMetadata(meth_gr))), nrow = 1))
      colnames(fix_metadata) = colnames(elementMetadata(meth_gr))
      elementMetadata(overhang_fix) = fix_metadata
      meth_gr = append(overhang_fix, meth_gr)
    }

    BottomStrandToFix = ifelse(context == "CG", "+", "-")
    if (as.character(strand(meth_gr[length(meth_gr)]))==BottomStrandToFix){
      print("Strand collapsing: Fixing bottom (right) overhang")
      overhang_fix = GRanges(seqnames = unique(as.character(seqnames(meth_gr))), ranges = IRanges(start(meth_gr[length(meth_gr)])+1, width = 1), strand = ifelse(context == "CG", "-", "+"))
      fix_metadata = data.frame(matrix(0, ncol = length(colnames(elementMetadata(meth_gr))), nrow = 1))
      colnames(fix_metadata) = colnames(elementMetadata(meth_gr))
      elementMetadata(overhang_fix) = fix_metadata
      meth_gr = append(meth_gr, overhang_fix)
    }

    meth_gr_collapsed <- meth_gr[seq(ifelse(context == "CG", 1, 2),length(meth_gr),by=2)]
    if (context == "CG"){
      end(meth_gr_collapsed) <- end(meth_gr_collapsed)+1
    } else if (context == "GC"){
      start(meth_gr_collapsed) <- start(meth_gr_collapsed)-1
    }

    values(meth_gr_collapsed) <- as.matrix(values(meth_gr[seq(1,length(meth_gr),by=2)]))+as.matrix(values(meth_gr[seq(2,length(meth_gr),by=2)]))

  } else {

    meth_gr_collapsed = GRanges()

  }

  return(meth_gr_collapsed)

}

CoverageFilter <- function(meth_gr, thr, context){

  Tcounts = grep('_T\\>',colnames(elementMetadata(meth_gr)))
  Mcounts = grep('_M\\>',colnames(elementMetadata(meth_gr)))

  MethylationMatrix = as.matrix(elementMetadata(meth_gr)[,Mcounts])/as.matrix(elementMetadata(meth_gr)[,Tcounts])
  #filter for coverage
  CovFilter=as.matrix(elementMetadata(meth_gr)[,Tcounts])>thr
  for (i in seq_len(ncol(MethylationMatrix))){
    MethylationMatrix[!CovFilter[,i],i] = NA
    }
  #bind the GRanges with the scores
  elementMetadata(meth_gr) = NULL
  meth_gr_filtered = resize(meth_gr, 1, fix=ifelse(context == "CG", 'start', "end"))
  elementMetadata(meth_gr_filtered)$MethylationFraction = MethylationMatrix

  return(meth_gr_filtered)

}

#' Call Context Methylation
#' @export
#'
#' @param meth_gr Methylation GRange object as returned by QuasR function qMeth
#' @param c0 coverage threshold
#' @param genome Primary sequence of genome of interest. E.g. for BSgenome.Mmusculus.UCSC.mm10 this would be Mmusculus
CallContextMethylation=function(meth_gr,cO,genome){

  # Subset Cytosines by genomic context
  meth_gr_CGs = findContextCytosines(meth_gr, genome, "CG")
  meth_gr_GCs = findContextCytosines(meth_gr, genome, "GC")

  # collapse strands
  meth_gr_CGs_collapsed = CollapeStrands(meth_gr_CGs, "CG")
  meth_gr_GCs_collapsed = CollapeStrands(meth_gr_GCs, "GC")

  #filter for coverage
  meth_gr_CGs_collapsed_filtered = CoverageFilter(meth_gr_CGs_collapsed, thr = cO, context = "CG")
  meth_gr_GCs_collapsed_filtered = CoverageFilter(meth_gr_GCs_collapsed, thr = cO, context = "GC")

  #disclose context
  GCGs=as.matrix(findOverlaps(meth_gr_GCs_collapsed_filtered,meth_gr_CGs_collapsed_filtered,type='equal'))
  meth_gr_CGs_collapsed_filtered$type='CGH'
  meth_gr_GCs_collapsed_filtered$type='GCH'
  meth_gr_CGs_collapsed_filtered$type[GCGs[,2]]='GCG'
  meth_gr_GCs_collapsed_filtered$type[GCGs[,1]]='GCG'

  meth_gr_combined = list(CG = meth_gr_CGs_collapsed_filtered, GC = meth_gr_GCs_collapsed_filtered)

  return(meth_gr_combined)

  }


