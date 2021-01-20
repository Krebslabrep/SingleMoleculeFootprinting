string.split=function(string,sep,pos){unlist(lapply(string,function(x){lapply(strsplit(x,sep),'[',pos)}))}



# used to be in single molecule manipulation functions
VectorizeReads=function(target,mergedMat){ #vectorize reads from matrices in order to plot in base space

  GCmet=as.vector(t(mergedMat))
  GCstartP=rep(as.numeric(colnames(mergedMat)),nrow(mergedMat))
  GCreadID=unlist(lapply(seq(nrow(mergedMat)),function(i){rep(i,ncol(mergedMat))}))
  GCgr=GRanges(seqnames(target),IRanges(as.numeric(colnames(mergedMat)),as.numeric(colnames(mergedMat))))
  list(GCstartP,GCreadID,GCmet)
}

#' Find genomic Cytosines for given context across given chromosome
#'
#' @param context Context of interest
#' @param chr chromsome of interest
#'
#' @return ??
#'
findCytosinesInContext = function(context, chr){
  cytosines <- matchPattern(context,Mmusculus[[chr]], fixed = ifelse(context == "DGCHN", F, T))
  cytosines_range <- GRanges(seqnames = Rle(chr),ranges = IRanges(start(cytosines), end = end(cytosines)),strand="+")
  if(context == "DGCHN"){
    cytosines_range <-resize(cytosines_range,2,fix='center')
  }
  if(context == "GC" || context == "DGCHN"){
    start(cytosines_range) <- start(cytosines_range)+1
  } else if (context == "CG"){
    end(cytosines_range) <- end(cytosines_range)-1
  }
  start(cytosines_range)
}

CollapseStrands_SelectCytosines = function(methAln, regDF){

  if(length(grep("_NO_",regDF[1,4])) > 0){
    #GC
    print("Detected experiment type: NO")
    CytosinesInContext = findCytosinesInContext(context = "DGCHN", chr = as.character(regDF$seqnames))
    methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]+1
    selCs <- methAln$Cid %in% CytosinesInContext
  }else if(length(grep("_SS_",regDF[1,4])) > 0){
    #CG
    print("Detected experiment type: SS")
    CytosinesInContext = findCytosinesInContext(context = "CG", chr = as.character(regDF$seqnames))
    methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]-1
    selCs <- methAln$Cid %in% CytosinesInContext
  }else if(length(grep("_DE_",regDF[1,4])) > 0){
    #GC&CG simultaneou1y. strand collapse is more complicated. prioritize GC in GCG minus strand collapse
    print("Detected experiment type: DE")
    CytosinesIn_CG = findCytosinesInContext(context = "CG", chr = as.character(regDF$seqnames))
    CytosinesIn_GC = findCytosinesInContext(context = "GC", chr = as.character(regDF$seqnames))
    plusStrandSel <- methAln$strand=="+"
    SelOvGC_plus <- ((methAln$Cid %in% CytosinesIn_GC) & plusStrandSel)
    SelOvGC_minus <- ((methAln$Cid+1) %in% CytosinesIn_GC) & !plusStrandSel
    SelOvCG_plus <- ((methAln$Cid %in% CytosinesIn_CG) & plusStrandSel)
    SelOvCG_minus <- ((methAln$Cid-1) %in% CytosinesIn_CG) & !plusStrandSel
    methAln$Cid[SelOvGC_minus] <- methAln$Cid[SelOvGC_minus]+1
    methAln$Cid[!SelOvGC_minus & SelOvCG_minus] <- methAln$Cid[!SelOvGC_minus & SelOvCG_minus]-1
    selCs <- SelOvGC_plus | SelOvGC_minus | SelOvCG_plus | SelOvCG_minus
  }else{stop("Sample name error. No _McvPI_ , _SsssI_ or _McvPISsssI_")}

  return(selCs)

}

Filter_ConversionRate_and_Duplicates = function(methAln, reg, conv.rate, remove.duplicates){

  #WCW for conversion and de-duplication
  dm3_WCWs_chr_XSV_p <- matchPattern("HCH",Mmusculus[[as.character(seqnames(reg))]],fixed="subject")
  dm3_WCWs_chr_XSV_n <- matchPattern(reverseComplement(DNAString("HCH")),Mmusculus[[as.character(seqnames(reg))]],fixed="subject")

  dm3_WCWs_chr_p <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_WCWs_chr_XSV_p), end = end(dm3_WCWs_chr_XSV_p)),strand="+")
  dm3_WCWs_chr_n <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_WCWs_chr_XSV_n), end = end(dm3_WCWs_chr_XSV_n)),strand="-")

  #start(dm3_WCWs_chr) <- start(dm3_WCWs_chr)+1
  dm3_WCWs_chr_coord <- c(start(dm3_WCWs_chr_p)+1,start(dm3_WCWs_chr_n)+1)

  ###################
  #identify molecules with low conversion rates
  ###################
  methAln_WCW<-methAln
  selCs_WCW <- methAln_WCW$Cid %in% dm3_WCWs_chr_coord
  readIDs=methAln_WCW[[1]][selCs_WCW]
  metV=methAln_WCW[[4]][selCs_WCW]
  molConvRates=tapply(metV,readIDs,function(x){
    mean(x,na.rm=T)
  })
  conversion.id=(1-molConvRates)*100>conv.rate
  #
  ###################
  #remove duplicates based on conversion errors (useful for amplicon data)
  ###################
  #identify duplicates
  if(remove.duplicates==T){
    sm=split(metV,readIDs)
    sm2=lapply(seq_along(sm),function(i){paste(sm[[i]],sep='',collapse ='')})
    unique.id=!duplicated(sm2)
    selReads=names(conversion.id)[conversion.id&unique.id]
    methAlnF=lapply(seq_along(methAln),function(i){methAln[[i]][methAln[['aid']]%in% selReads]})

  }else{
    selReads=names(conversion.id)[conversion.id]
    methAlnF=lapply(seq_along(methAln),function(i){methAln[[i]][methAln[['aid']]%in% selReads]})
  }

}

BinMethylation = function(Crange, st, binCoord){

  midP = start(resize(st,1,fix='center'))
  binRange = GRanges(seqnames(st), IRanges(ifelse(strand(st) =='+', midP+binCoord[1], midP-binCoord[2]),
                                           ifelse(strand(st)=='+', midP+binCoord[2], midP-binCoord[1])))
  binOverlaps = findOverlaps(ranges(binRange), Crange)
  #compute the methylation vectors
  binMethylation = elementMetadata(Crange)$meth[subjectHits(binOverlaps)]
  binReads = paste(names(st)[queryHits(binOverlaps)], elementMetadata(Crange)$readID[subjectHits(binOverlaps)], sep='_')
  #group Cs per region/per_read
  grouped_binMethylation = round(tapply(binMethylation,binReads,mean))
  unique_binReads = sort(unique(binReads))

  return(list(grouped_binMethylation, unique_binReads))

}

# used to be in plotting
OneTFstates = function(){

  allPos=expand.grid(c(0,1),c(0,1),c(0,1))
  patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  #using only 'pure' states
  states=list(
    bound=patternStrings[6],
    accessible=patternStrings[8],
    closed=patternStrings[1],
    unassigned=patternStrings[!seq_along(patternStrings)%in%c(1,6,8)]
  )
  states = rev(states)

  return(states)

}

TFpairStates = function(){
  #define states for each factor separately
  allPos=expand.grid(c(0,1),c(0,1))
  TF1=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  names(TF1)=c('nucleosome','unassigned','bound','accessible')
  TF2=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  names(TF2)=c('nucleosome','bound','unassigned','accessible')

  #create a combined state factor
  combined_states=unlist(lapply(seq_along(TF1),function(i){
    lapply(seq_along(TF2),function(j){
      paste(TF1[i],TF2[j],sep='')
    })}))
  names(combined_states)=unlist(lapply(seq_along(TF1),function(i){
    lapply(seq_along(TF2),function(j){
      paste(names(TF1)[i],names(TF2)[j],sep='_')
    })}))
  combined_statesF=as.factor(unlist(lapply(seq_along(combined_states),function(i){rep(names(combined_states[i]),length(combined_states[[i]]))}))[order(unlist(combined_states))])
  combined_statesF=factor(combined_statesF,levels=names(combined_states))
  stateM=cbind(string.split(as.character(combined_statesF),'_',1),string.split(as.character(combined_statesF),'_',2))
  grouped_states=list(
    combined_states[stateM[,1]=='bound'& stateM[,2]=='bound'],
    combined_states[stateM[,1]=='bound'& !stateM[,2]=='bound' & !stateM[,2]=='nucleosome'],
    combined_states[!stateM[,1]=='bound'& !stateM[,1]=='nucleosome' &stateM[,2]=='bound'],
    combined_states[(!stateM[,1]=='bound'& !stateM[,2]=='bound')|(stateM[,1]=='bound'& stateM[,2]=='nucleosome')|(stateM[,2]=='bound'& stateM[,1]=='nucleosome')])
  grouped_states = rev(grouped_states)

  return(grouped_states)
}

SortReads_internal = function(SortedReads, SM_mat, isClusters){

  read_sort=SortedReads[[1]]
  read_sort=read_sort[!is.na(names(read_sort))]
  if(isClusters){
    states = TFpairStates()
    stN=string.split(names(read_sort),'_',4)
  } else {
    states = OneTFstates()
    stN=unlist(lapply(names(read_sort), function(x){strsplit(x, "_")[[1]][3]}))
  }
  names(read_sort)=stN
  read_sort.e=vector("list", length(unlist(states)))
  names(read_sort.e)=unlist(states)
  read_sort.e[stN]=read_sort
  read_sort=read_sort.e
  read_sort=lapply(seq_along(read_sort),function(i){read_sort[[i]][read_sort[[i]]%in%rownames(SM_mat)]}) # I don't get it, these should already be the same
  names(read_sort)=unlist(states)
  read_sort=read_sort[(unlist(states))]
  read_sort
  return(read_sort)

}
