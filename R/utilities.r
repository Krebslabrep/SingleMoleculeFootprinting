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
    closed=patternStrings[1],
    unassigned=patternStrings[!seq_along(patternStrings)%in%c(1,6,8)]
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

  return(grouped_states)
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

