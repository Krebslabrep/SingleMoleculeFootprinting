string.split=function(string,sep,pos){unlist(lapply(string,function(x){lapply(strsplit(x,sep),'[',pos)}))}

# Useful for plotting
VectorizeReads=function(target,mergedMat){ #vectorize reads from matrices in order to plot in base space

  GCmet=as.vector(t(mergedMat))
  GCstartP=rep(as.numeric(colnames(mergedMat)),nrow(mergedMat))
  GCreadID=unlist(lapply(seq(nrow(mergedMat)),function(i){rep(i,ncol(mergedMat))}))
  GCgr=GRanges(seqnames(target),IRanges(as.numeric(colnames(mergedMat)),as.numeric(colnames(mergedMat))))
  list(GCstartP,GCreadID,GCmet)
}

#' Design states for single TF case
#'
#' @return list of states
#'
OneTFstates = function(){

  allPos=expand.grid(c(0,1),c(0,1),c(0,1))
  patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  #using only 'pure' states
  states=list(
    bound=patternStrings[6],
    accessible=patternStrings[8],
    closed=patternStrings[c(1,2,5)],
    unassigned=patternStrings[!seq_along(patternStrings)%in%c(1,2,5,6,8)]
  )

  return(states)

}

#' Design states for TF pair case
#'
#' @return list of states
#'
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
  names(grouped_states) = c("misc", "misc_bound", "bound_misc", "bound_bound")

  return(grouped_states)
}

#' Design states for promoters
#'
#' @return list of states
#'
Promoterstates = function(){
  
  allPos=expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
  patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  #using only 'pure' states
  states=list(
    unassigned=patternStrings[!seq_along(patternStrings) %in% c(1:5,8,9,10,12,13,14:16)],
    nucleosome=patternStrings[c(1:5,9,13)],
    unbound=patternStrings[c(8,15,16)],
    PIC=patternStrings[12],
    PIC.polII=patternStrings[10],
    polII=patternStrings[14]
  )
  
  return(states)
}


SortReads_internal = function(SortedReads, SM_mat, isClusters){ # this orders readIDs based on arbotrary states order

  read_sort=SortedReads#[[1]]
  read_sort=read_sort[!is.na(names(read_sort))]
  if(isClusters){
    states = TFpairStates()
    stN=string.split(names(read_sort),'_',4)
  } else {
    states = OneTFstates()
    stN=names(read_sort)#unlist(lapply(names(read_sort), function(x){strsplit(x, "_")[[1]][3]}))
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

# Utility functions to SampleCorrelation function
panel.jet <- function(...) {
  jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  smoothScatter(..., nrpoints=0, add=TRUE, colramp=jet.colors) }

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
}

#' @importFrom IRanges cor
#' @importFrom stats cor.test symnum
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- (cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

  test <- cor.test(x,y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))

  text(0.5, 0.5, txt, cex=0.6*cex)

}

#' Calculate colMeans after dropping zeros
#' 
#' @param MethSM one single molecule sparse matrix
#' 
#' @import Matrix
#'  
#' @return colMeans (N.b. this is +1 based)
#'  
colMeans_drop0 <- function (MethSM) {
  nnz_per_col <- diff(MethSM@p)
  # nnz_per_col[nnz_per_col == 0] <- 1  ## just avoid doing 0 / 0
  return(Matrix::colSums(MethSM) / nnz_per_col)
}

#' Calculate rowMeans after dropping zeros
#' 
#' @param MethSM one single molecule sparse matrix
#' 
#' @import Matrix
#'  
#' @return rowMeans (N.b. this is +1 based)
#'  
rowMeans_drop0 <- function (MethSM) {
  RowInd <- MethSM@i + 1
  nnz_per_row <- tabulate(RowInd)
  if (length(nnz_per_row) < MethSM@Dim[1]) {
    nnz_per_row = c(nnz_per_row, rep(0, MethSM@Dim[1] - length(nnz_per_row)))
  }
  # nnz_per_row[nnz_per_row == 0] <- 1  ## just avoid doing 0 / 0
  return(Matrix::rowSums(MethSM) / nnz_per_row)
}

#' Recalculate *_T and *_M values in MethGR object after filtering reads e.g. for conversion rate
#' 
#' @param MethGR GRanges object of methylation call
#' @param MethSM Single Molecule methylation matrix
#' @param MethSM_filtered Single Molecule methylation matrix after filtering reads
#' @param sampleIndex index for sample to treat. It serves as a correspondence between the index of the SM matrix and the order samples appear in the elementMetadata() columns
#' 
#' @import Matrix
#' @import GenomicRanges
#' 
#' @return MethGR with recalculated counts
#' 
filter_reads_from_MethGR = function(MethGR, MethSM, MethSM_filtered, sampleIndex){
  
  DiscardedReads = MethSM@Dimnames[[1]][!(MethSM@Dimnames[[1]] %in% MethSM_filtered@Dimnames[[1]])]
  DiscardedReads = ifelse(length(DiscardedReads) == 0, 0, DiscardedReads)
  AffectedCytosines = MethSM[DiscardedReads,,drop=FALSE]@Dimnames[[2]][colSums(MethSM[DiscardedReads,,drop=FALSE]) > 0]
  AffectedCytosines = ifelse(length(AffectedCytosines) == 0, 0, AffectedCytosines)
  T_counts = diff(MethSM_filtered[,AffectedCytosines,drop=FALSE]@p)
  M_counts = colMeans_drop0(MethSM_filtered[,AffectedCytosines,drop=FALSE]) - 1
  elementMetadata(MethGR)[start(MethGR) %in% AffectedCytosines, grep("_T$", colnames(elementMetadata(MethGR)), value=TRUE)[sampleIndex]] = T_counts
  elementMetadata(MethGR)[start(MethGR) %in% AffectedCytosines, grep("_M$", colnames(elementMetadata(MethGR)), value=TRUE)[sampleIndex]] = M_counts
  return(MethGR)
  
}

#' Implementation performing a similar operation of plyr::rbind.fill.matrix but for sparseMatrix
#' 
#' @param x sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' @param y sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' 
#' @details N.b. only possible fill at the moment is 0
#' 
#' @importFrom Matrix rsparsematrix
#' 
#' @export 
#' 
rbind.fill.matrix.sparse = function(x,y){
  
  ymiss = colnames(x)[which(is.na(match(colnames(x),colnames(y))))]
  ybind = Matrix::rsparsematrix(nrow=as.double(nrow(y)),ncol=as.double(length(ymiss)),density = 0)
  colnames(ybind)<-ymiss
  
  xmiss = colnames(y)[which(is.na(match(colnames(y),colnames(x))))]
  xbind = Matrix::rsparsematrix(nrow=as.double(nrow(x)),ncol=as.double(length(xmiss)),density = 0)
  colnames(xbind) = xmiss

  if (ncol(xbind)>0){
    x = cbind2(x,xbind)
    x = x[,order(colnames(x)),drop=FALSE]
  }
  if(ncol(ybind)>0){
    y = cbind2(y,ybind)
    y = y[,order(colnames(y)),drop=FALSE]
  }
  
  result = rbind2(x,y[,order(match(colnames(y),colnames(x))),drop=FALSE])
  if (all(result@Dim > 0)){
    rownames(result) = c(rownames(x), rownames(y)) # for some reason rbind2 drop rownames
    # result = result[,order(colnames(result), decreasing = FALSE)] # This shouldn't be necessary for rows
  }
  
  return(result)
  
}

#' Implementation performing a similar operation of rbind.fill.matrix.sparse but for columns
#' 
#' @param x sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' @param y sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' 
#' @details N.b. only possible fill at the moment is 0
#' 
#' @export 
#' 
cbind.fill.matrix.sparse = function(x,y){
  
  ymiss = rownames(x)[which(is.na(match(rownames(x),rownames(y))))]
  ybind = rsparsematrix(nrow=as.double(length(ymiss)),ncol=as.double(ncol(y)),density = 0)
  rownames(ybind) = ymiss
  
  xmiss = rownames(y)[which(is.na(match(rownames(y),rownames(x))))]
  xbind = rsparsematrix(nrow=as.double(length(xmiss)),ncol=as.double(ncol(x)),density = 0)
  rownames(xbind) = xmiss
  
  x = rbind2(x,xbind)
  y = rbind2(y,ybind)
  
  result = cbind2(x,y[order(match(rownames(y),rownames(x))),,drop=FALSE])
  if (all(result@Dim > 0)){
    colnames(result) = c(colnames(x), colnames(y)) # for some reason cbind2 drop colnames
    result = result[,order(colnames(result), decreasing = FALSE),drop=FALSE]
  }
  
  return(result)
  
}
