#' Utility to getGCMatrix
#'
fuseReadMat <- function(resp,resm){ #assumes that reads are no on different strands
  uReads=unique(c(rownames(resp),rownames(resm)))
  uCs=unique(c(colnames(resp),colnames(resm)))
  matCGo=matrix(nrow=length(uReads),ncol=length(uCs))
  colnames(matCGo)=uCs
  rownames(matCGo)=uReads
  matCGo[rownames(resp),colnames(resp)]=resp #fill up positive reads
  #add negative reads
  NAi=apply(resm,1,function(x){sum(is.na(x))==length(x)})
  matCGo[rownames(resm[!NAi,]),colnames(resm[!NAi,])]=resm[!NAi,]
  matCGo
}

#' Get Single Molecule methylation matrix
#'
#' Used internally as the first step in getGCMatrix
#'
#' @param QuasRprj QuasR project object as returned by calling the QuasR function qAling on previously aligned data
#' @param range GenimocRange representing the genomic region of interest
#' @param sample One of the sample names as reported in the SampleName field of the QuasR pointer file provided to qAlign. N.b. all the files with the passed sample name will be used to call methylation
#'
#' @return Single molecule methylation matrix (all Cytosines)
#'
#' @export
getCMethMatrix<-function(QuasRprj,range,sample){

  Cs=QuasR::qMeth(QuasRprj[grep(sample, QuasRprj@alignments$SampleName)], query=range, mode="allC",reportLevel="alignment", collapseBySample = T)
  # use data.table to get a 1,0 matrix of methylation profiles
  # all.cids=unique(Cs[[sample]]$Cid) # get all possible C locations
  # make the data.table object
  dt=data.table::data.table(meth=Cs[[sample]]$meth ,aid=Cs[[sample]]$aid ,cid=Cs[[sample]]$Cid)
  spread_dt = data.table::dcast(dt, aid ~ cid, value.var = "meth")
  meth_matrix = as.matrix(spread_dt, rownames = "aid")
  return(meth_matrix)

}

#' Extract context cytosines from SM methylation matrix
#'
#' Supposed to be used after getCMethMatrix.
#'
#' @param QuasRprj QuasR project object as returned by calling the QuasR function qAling on previously aligned data
#' @param range GenimocRange representing the genomic region of interest
#' @param sample One of the sample names as reported in the SampleName field of the QuasR pointer file provided to qAlign. N.b. all the files with the passed sample name will be used to call methylation
#' @param genome BS genome
#' @param conv.rate default = 80
#' @param destrand defualt = TRUE
#'
#' @return Single molecule methylation matrixes (contexts resolved)
#'
#' @export
getGCMatrix<-function(QuasRprj, range, sample, genome, conv.rate=80, destrand=TRUE){

  matC = getCMethMatrix(QuasRprj, range, sample)
  chr = seqnames(range)

  Cpos=as.numeric(colnames(matC))
  rGR=GRanges(rep(chr,length(Cpos)),IRanges(Cpos,Cpos))
  rSeq=getSeq(genome,resize(rGR,3,fix='center'))
  GCpos=vcountPattern('GC',rSeq)==1
  CGpos=vcountPattern('CG',rSeq)==1
  #       GCGpos=vcountPattern('GCG',rSeq)
  GCi= GCpos #& !GCGpos
  CGi= CGpos #& !GCGpos
  Ci= Cpos & !GCpos & !CGpos #& !GCGpos

  ########
  #Filter bases
  ########

  # filter based on conversion
  ConvRate=100-(rowMeans(matC[,Ci],na.rm=T)*100)
  Convi=ConvRate>conv.rate
  matGC=matC[Convi,GCi,drop=FALSE] # get matrices
  matCG=matC[Convi,CGi,drop=FALSE]
  ########
  #Destrand
  ########
  # bases G or C
  # will help determine the strands
  bpsGC=as.character(getSeq(genome,rGR))[GCi]
  bpsCG=as.character(getSeq(genome,rGR))[CGi]

  if(destrand){
    #Destrand GC
    # column names that are for plus and minus strand
    colMinus=(colnames(matGC))[bpsGC=="G"]
    colPlus= (colnames(matGC))[bpsGC=="C"]

    # get the minus strand columns
    # remove empty rows
    # order by column names
    # add +1 to column names
    matGCM=matGC[, colnames(matGC) %in% colMinus,drop=FALSE]
    matGCM=matGCM[rowSums(!is.na(matGC))>0,,drop=FALSE]
    matGCM=matGCM[,order(as.numeric(colnames(matGCM))),drop=FALSE]
    colnames(matGCM)=as.numeric(colnames(matGCM))+1

    # get the plus strand columns
    # remove empty rows
    # order by column names
    matGCP=matGC[ ,colnames(matGC) %in% colPlus,drop=FALSE]
    matGCP=matGCP[rowSums(!is.na(matGCP))>0,,drop=FALSE]
    matGCP=matGCP[,order(as.numeric(colnames(matGCP))),drop=FALSE] # reorder columns

    # temp matrices for plus and minus strand reads
    resp=matrix(NA,ncol=length(unique(c(colnames(matGCM),colnames(matGCP)))),
                nrow=nrow(matGCP))
    rownames(resp)=rownames(matGCP)
    resm=matrix(NA,ncol=length(unique(c(colnames(matGCM),colnames(matGCP)))),
                nrow=nrow(matGCM))
    rownames(resm)=rownames(matGCM)
    # get col names as base numbers
    cols=unique(c(colnames(matGCM),colnames(matGCP)))
    cols=cols[order(as.numeric(cols))]
    colnames(resp)=cols
    colnames(resm)=cols

    # populate temporary matrices
    resp[,colnames(resp) %in% colnames(matGCP)]=matGCP
    resm[,colnames(resm) %in% colnames(matGCM)]=matGCM
    rownames(resp)=rownames(matGCP)
    rownames(resm)=rownames(matGCM)
    #add the read names


    matGC=fuseReadMat(resp,resm)#cbind(resp,resm)

    matGC=matGC[,order(as.numeric(colnames(matGC))),drop=FALSE] # reorder columns

    #Destrand CG
    # column names that are for plus and minus strand
    colMinus=(colnames(matCG))[bpsCG=="G"]
    colPlus= (colnames(matCG))[bpsCG=="C"]

    # get the minus strand columns
    # remove empty rows
    # order by column names
    # add -1 to column names
    matCGM=matCG[, colnames(matCG) %in% colMinus,drop=FALSE]
    matCGM=matCGM[rowSums(!is.na(matCG))>0,,drop=FALSE]
    matCGM=matCGM[,order(as.numeric(colnames(matCGM))),drop=FALSE]
    colnames(matCGM)=as.numeric(colnames(matCGM))-1

    # get the plus strand columns
    # remove empty rows
    # order by column names
    matCGP=matCG[ ,colnames(matCG) %in% colPlus,drop=FALSE]
    matCGP=matCGP[rowSums(!is.na(matCGP))>0,,drop=FALSE]
    matCGP=matCGP[,order(as.numeric(colnames(matCGP))),drop=FALSE] # reorder columns

    # temp matrices for plus and minus strand reads
    resp=matrix(NA,ncol=length(unique(c(colnames(matCGM),colnames(matCGP)))),nrow=nrow(matCGP))
    resm=matrix(NA,ncol=length(unique(c(colnames(matCGM),colnames(matCGP)))),nrow=nrow(matCGM))
    rownames(resp)=rownames(matCGP)
    rownames(resm)=rownames(matCGM)
    # get col names as base numbers
    cols=unique(c(colnames(matCGM),colnames(matCGP)))
    cols=cols[order(as.numeric(cols))]
    colnames(resp)=cols
    colnames(resm)=cols

    # populate temporary matrices
    resp[,colnames(resp) %in% colnames(matCGP)]=matCGP
    resm[,colnames(resm) %in% colnames(matCGM)]=matCGM
    rownames(resp)=rownames(matCGP)
    rownames(resm)=rownames(matCGM)
    #merge the two matrices

    matCG=fuseReadMat(resp,resm)

    matCG=matCG[,order(as.numeric(colnames(matCG))),drop=FALSE] # reorder columns

  }else{

    matGC=matGC[,order(as.numeric(colnames(matGC))),drop=FALSE] # reorder columns
    matCG=matCG[,order(as.numeric(colnames(matCG))),drop=FALSE]
  }
  ########
  #Remove the emty row and collumns
  ########


  matGC=matGC[,colSums(!is.na(matGC))>0,drop=FALSE] # remove no coverage columns
  matCG=matCG[,colSums(!is.na(matCG))>0,drop=FALSE]
  matGC=matGC[rowSums(!is.na(matGC))>0,,drop=FALSE] # remove no coverage rows
  matCG=matCG[rowSums(!is.na(matCG))>0,,drop=FALSE]
  list(matGC=matGC,matCG=matCG)
}

#' Merge context resolved SM matrixes
#'
#' Supposed to be used after getGCMatrix
#'
#' @param GC_CG_mat Output of getGCMatrix
#'
#' @return Merged single molecule methylation matrix (contexts resolved)
#'
#' @export
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

#' Filter Cytosines in context
#'
#' @param meth_gr Granges obj of average methylation
#' @param genome BSgenome
#' @param context Context of interest (e.g. "GC", "CG",..)
#'
#' @return filtered Granges obj
#'
FilterContextCytosines <- function(meth_gr, genome, context){

  chr = as.character(unique(seqnames(meth_gr)))
  GenomicCytosines = matchPattern(context, genome[[chr]])
  GenomicCytosines_gr = GRanges(seqnames = chr, ranges = GenomicCytosines)
  meth_gr[meth_gr %over% GenomicCytosines_gr]

}

#' Collapse strands
#'
#' @param meth_gr
#' @param context
#'
#' @return meth_gr with collapsed strands (everything turned to - strand)
#'
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

#' Filter Cs for coverage
#'
#' @param meth_gr
#' @param thr
#' @param context
#'
#' @return filtered meth_gr
#'
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
#'
#' Can deal with multiple samples
#'
#' @export
#'
#' @param QuasRprj QuasR project object as returned by calling the QuasR function qAling on previously aligned data
#' @param range GenimocRange representing the genomic region of interest
#' @param coverage coverage threshold. Defaults to 20.
#' @param genome Primary sequence of genome of interest. E.g. for BSgenome.Mmusculus.UCSC.mm10 this would be Mmusculus
#'
#' @return List with two Granges objects: GC and CG average methylation
#'
CallContextMethylation=function(QuasRprj, range, coverage=20, genome){

  meth_gr <- QuasR::qMeth(QuasRprj, mode="allC", range, collapseBySample = T)

  # Subset Cytosines by genomic context
  meth_gr_CGs = FilterContextCytosines(meth_gr, genome, "CG")
  meth_gr_GCs = FilterContextCytosines(meth_gr, genome, "GC")

  # collapse strands
  meth_gr_CGs_collapsed = CollapeStrands(meth_gr_CGs, "CG")
  meth_gr_GCs_collapsed = CollapeStrands(meth_gr_GCs, "GC")

  #filter for coverage
  meth_gr_CGs_collapsed_filtered = CoverageFilter(meth_gr_CGs_collapsed, thr = coverage, context = "CG")
  meth_gr_GCs_collapsed_filtered = CoverageFilter(meth_gr_GCs_collapsed, thr = coverage, context = "GC")

  #disclose context
  GCGs=as.matrix(findOverlaps(meth_gr_GCs_collapsed_filtered,meth_gr_CGs_collapsed_filtered,type='equal'))
  meth_gr_CGs_collapsed_filtered$type='CGH'
  meth_gr_GCs_collapsed_filtered$type='GCH'
  meth_gr_CGs_collapsed_filtered$type[GCGs[,2]]='GCG'
  meth_gr_GCs_collapsed_filtered$type[GCGs[,1]]='GCG'

  meth_gr_combined = list(CG = meth_gr_CGs_collapsed_filtered, GC = meth_gr_GCs_collapsed_filtered)

  return(meth_gr_combined)

  }


#' Single molecule methylation matrix
#'
#' needs to be run on each sample separately
#'
#'
#'
#' @export
SingleMoleculeMatrix = function(QuasRprj, st, sample, coverage=10, genome){

  print("Producing Single Molecule matrix")

  sm.mat=getCMethMatrix(QuasRprj,st,sample)

  print(sm.mat)

  if (dim(sm.mat)[2]>0 & dim(sm.mat)[1]>coverage){

    GC_CGmat=getGCMatrix(matC = sm.mat,chr=as.character(seqnames(st)[1]),genome=genome,destrand=TRUE, conv.rate = 80)
    #separate context
    if(length(grep("_NO_",regDF[1,4])) > 0){
      #GC
      print("Detected experiment type: NO")
      st.seq=getSeq(genome,st)
      GCpos=start(st)+start(matchPattern(DNAString("DGCHN"),st.seq[[1]],fixed="subject"))+1
      SM=GC_CGmat[['matGC']][,colnames(GC_CGmat[['matGC']])%in% as.character(GCpos)]

    }else if(length(grep("_SS_",regDF[1,4])) > 0){
      #CG
      print("Detected experiment type: SS")
      CGpos=start(st)+start(matchPattern(DNAString("NWCGW"),st.seq[[1]],fixed="subject"))+1
      SM=GC_CGmat[['matCG']][,colnames(GC_CGmat[['matCG']])%in% as.character(CGpos)]

    }else if(length(grep("_DE_",regDF[1,4])) > 0){
      #GC&CG simultaneously. strand collapse is more complicated. prioritize GC in GCG minus strand collapse
      print("Detected experiment type: DE")
      SM = cbind(GC_CGmat$matGC, GC_CGmat$matCG)

    }else{stop("Sample name error. No _McvPI_ , _SsssI_ or _McvPISsssI_")}

    return(SM)


  } else {stop(paste0("length(sm.mat)<0 || length(unlist(SortedReads))<", coverage))}

}
