#' Conversion rate
#'
#' calculate sequencing library conversion rate on a chromosome of choice
#' @param sampleSheet QuasR sample sheet
#' @param genome BS genome
#' @param chr chromosome to calculate conversion rate on (default: 19)
#'
#' @export
ConversionRate = function(sampleSheet, genome, chr=19){

  prj <- qAlign(sampleFile=sampleSheet,
                genome=genome@pkgname,
                aligner = "Rbowtie",
                paired="fr",
                bisulfite="undir",
                projectName="prj",
                checkOnly=F)
  prj@aligner = "Rbowtie"
  nr.cores = length(unique(prj@alignments$SampleName))

  # check convertion rate:
  # how many of the non CG/GC cytosines are methylated?

  seq_length = seqlengths(genome)
  chr = tileGenome(seq_length[chr], tilewidth=max(seq_length[chr]), cut.last.tile.in.chrom=TRUE)
  cl = makeCluster(nr.cores)
  methylation_calls_C = qMeth(prj, query = chr, mode="allC", reportLevel="C", keepZero = T, clObj = cl, asGRanges = T, collapseBySample = F)

  seqContext = BSgenome::getSeq(genome,resize(methylation_calls_C,3, fix='center'))
  GCc = Biostrings::vcountPattern(DNAString("GCN"),seqContext,fixed=F) # take all context for exclusion
  CGc = Biostrings::vcountPattern(DNAString("NCG"),seqContext,fixed=F) # only valid for calling conversion rates
  #SMF
  #non CG non GC
  out_c_met = methylation_calls_C[GCc==0 & CGc==0,]
  # CG_met = methylation_calls_C[!CGc==0,]
  # GC_met = methylation_calls_C[!GCc==0,]
  #non GC context
  tot.col = grep('R\\d_T',colnames(values(out_c_met)))
  tot.col = grep('_T$',colnames(values(out_c_met)))
  met.col = grep('_M$',colnames(values(out_c_met)))
  tot.c = as.matrix(values(out_c_met)[,tot.col])
  met.c = as.matrix(values(out_c_met)[,met.col])
  conv_rate = round((1-(colSums(met.c)/colSums(tot.c)))*100,1)

  stopCluster(cl)
  return(conv_rate)

}

#' Bait capture efficiency
#'
#' check bait capture efficiency. Expected to be ~70% for mouse genome
#' @param sampleSheet QuasR sample sheet
#' @param genome BS genome
#' @param baits Full path to bed file containing bait coordinates. If chromosome names are in e.g. "1" format, they'll be temporarily converted to "chr1"
#'
#' @export
BaitCapture = function(sampleSheet, genome, baits){

  prj <- QuasR::qAlign(sampleFile=sampleSheet,
                genome=genome@pkgname,
                aligner = "Rbowtie",
                paired="fr",
                bisulfite="undir",
                projectName="prj",
                checkOnly=F)
  prj@aligner = "Rbowtie"
  nr.cores = length(unique(prj@alignments$SampleName))

  if (length(grep("chr", seqlevels(BaitRegions))) == 0){
    seqlevels(BaitRegions) = paste0("chr", seqlevels(BaitRegions))
  }

  cl = makeCluster(nr.cores)
  InBaits=QuasR::qCount(prj,BaitRegions,clObj=cl)

  seq_length=seqlengths(genome)
  tiles <- tileGenome(seq_length, tilewidth=max(seq_length),cut.last.tile.in.chrom=TRUE)
  cl = makeCluster(nr.cores)
  GW=QuasR::qCount(prj,tiles,clObj=cl)

  capture_efficiency = c()
  for(n in 1:length(prj@alignments$SampleName)){
    capture_efficiency = c(capture_efficiency, sum(dplyr::select(dplyr::as_tibble(InBaits), -width)[,c(n)]) / sum(dplyr::select(dplyr::as_tibble(GW), -width)[,c(n)]))
  }

  stopCluster(cl)
  return(capture_efficiency)

}


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

# NEED TO DO: following function now doesn't work with external data, and correlation for example data is NA

#' Intersample correlation
#'
#' pair plot of sample correlations
#' @param sample Avg methylation object. Can also be set to "Example" to produce plot using example data of the kind specified by the @param CellType
#' @param context one of "AllCs", "DGCHN", "NWCGW". The first should be chosen for TKO experiments. For experiments carried on WT cells, we recommend checking both the "DGCHN" and "NWCGW" contexts by running this function once per context.
#' @param CellType Cell type to compare your samples to. At the moment, this can be one of "ES", "NP", "TKO".
#' @param saveAs Full path to output plot file
#'
#' @export
SampleCorrelation = function(samples, context, CellType, saveAs=NULL){

  # Get methylation data from previous SMF experiments
  contextmethylation <- readRDS("/g/krebs/krebs/analysis/SMF/MM/methCall/Context_methylation_call_20_all_samples.txt.rds")
  AllC=contextmethylation[[1]] # GRanges with all methylation sites
  metMat_ref=contextmethylation[[2]] # all methylation data
  remove(contextmethylation)

  # Pick reference data
  if (length(grep(CellType, colnames(metMat_ref))) > 0){
    refMeth = rowMeans(metMat_ref[,c(grep(CellType, colnames(metMat_ref))[1:2])])
  } else {
    print("Unsupported cell type, skipping comparison to reference dataset")
  }

  # Pick example data
  if (samples == "Example"){
    contextMet = metMat_ref[,c(grep(CellType, colnames(metMat_ref))[3])]
  } else {
    # prepare samples to be a matrix with the right coordinates
    # for (i in 1:ncol(contextMet)){
    #   ## Overlap the metCall object with the AllC to define its "absolute" position
    #   ov <- as.matrix(findOverlaps(contextMet[,i], AllC, type = "equal", ignore.strand =F))
    #   ##Use the overlap to fill in the correct positions in the matrix
    #   comparison_mat[,i + 1] [ov[,2]] <- contextMet$met_rate[ov[,1]]
    # }
  }

  # now we don't need this anymore
  remove(metMat_ref)

  # Prepare overall comparison matrix and add samples
  comparison_mat <- matrix(NA, ncol = ncol(contextMet)+1, nrow = length(refMeth))
  comparison_mat[,1] <- refMeth
  comparison_mat[comparison_mat=="NaN"] <- NA
  comparison_mat[,1:ncol(contextMet)] = contextMet
  colnames(comparison_mat) = c(CellType, colnames(contextMet))
  remove(contextMet)

  # Keep track of the rows with only NAs
  NAindex_comparison <- rowSums(is.na(comparison_mat))!=ncol(comparison_mat)

  # Take care of the context
  if (context == "AllCs") {
    CytosinesIndex = rep(TRUE, length(AllC))
  } else if (context == "DGCHN" | context == "NWCGW"){
    CytosinesIndex = AllC$type == context
  } else {
    stop("Unsupported context")
  }

  if (!is.null(saveAs)){
    pdf(saveAs)
    pairs(comparison_mat[CytosinesIndex & NAindex_comparison,], upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = colnames(comparison_mat))
    dev.off()
  } else {
    pairs(comparison_mat[CytosinesIndex & NAindex_comparison,], upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = colnames(comparison_mat))
  }

}
