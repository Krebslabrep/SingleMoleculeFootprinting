#' Conversion rate
#'
#' calculate sequencing library conversion rate on a chromosome of choice
#' @param sampleSheet QuasR sample sheet
#' @param genome BS genome
#' @param chr chromosome to calculate conversion rate on (default: 19)
#' @param cores number of cores for parallel processing. Defaults to 1
#'
#' @importFrom QuasR qMeth
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom parallel makeCluster stopCluster
#' @importFrom BSgenome getSeq
#' @importFrom IRanges resize
#' @importFrom Biostrings vcountPattern
#' @importFrom BiocGenerics grep
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#'
#' ConversionRatePrecision = ConversionRate(sampleSheet = Qinput, genome = BSgenome.Mmusculus.UCSC.mm10, chr = 19, cores = 1)
#'
ConversionRate = function(sampleSheet, genome, chr=19, cores=1){

  QuasRprj = GetQuasRprj(sampleSheet, genome)

  # check convertion rate:
  # how many of the non CG/GC cytosines are methylated?
  seq_length = seqlengths(genome)
  chr = tileGenome(seq_length[chr], tilewidth=max(seq_length[chr]), cut.last.tile.in.chrom=TRUE)
  cl = makeCluster(cores)
  methylation_calls_C = qMeth(QuasRprj, query = chr, mode="allC", reportLevel="C", keepZero = TRUE, clObj = cl, asGRanges = TRUE, collapseBySample = FALSE)
  stopCluster(cl)

  seqContext = getSeq(genome, trim(resize(methylation_calls_C, 3, fix='center')))
  GCc = vcountPattern(DNAString("GCN"), seqContext, fixed=FALSE) # take all context for exclusion
  CGc = vcountPattern(DNAString("NCG"), seqContext, fixed=FALSE) # only valid for calling conversion rates
  #SMF
  #non CG non GC
  out_c_met = methylation_calls_C[GCc==0 & CGc==0,]
  # CG_met = methylation_calls_C[!CGc==0,]
  # GC_met = methylation_calls_C[!GCc==0,]
  #non GC context
  # tot.col = grep('R\\d_T',colnames(values(out_c_met)))
  tot.col = grep('_T$',colnames(values(out_c_met)))
  met.col = grep('_M$',colnames(values(out_c_met)))
  tot.c = as.matrix(values(out_c_met)[,tot.col])
  met.c = as.matrix(values(out_c_met)[,met.col])
  conv_rate = round((1-(colSums(met.c)/colSums(tot.c)))*100,1)

  return(conv_rate)

}

#' Bait capture efficiency
#'
#' check bait capture efficiency. Expected to be ~70% for mouse genome
#' @param sampleSheet QuasR sample sheet
#' @param genome BS genome
#' @param baits GRanges obj of bait coordinates. We provide and example through SingleMoleculeFootprintingData::EnrichmentRegions_mm10.rds()
#' @param clObj cluster object to emply for parallel processing created using the parallel::makeCluster function. Defaults to NULL
#'
#' @import BiocGenerics
#' @importFrom QuasR qCount
#' @importFrom GenomeInfoDb seqlengths
#'
#' @return bait capture efficiency
#'
#' @export
#'
#' @examples
#' Qinput = paste0(tempdir(), "/NRF1Pair_Qinput.txt")
#' library(BSgenome.Mmusculus.UCSC.mm10)
#'
#' if(file.exists(Qinput)){
#'     # DO NOT RUN
#'     # clObj = parallel::makeCluster(5)
#'     # BaitRegions = SingleMoleculeFootprintingData::EnrichmentRegions_mm10.rds()
#'     # BaitCaptureEfficiency = BaitCapture(sampleSheet = Qinput, genome = BSgenome.Mmusculus.UCSC.mm10, baits = BaitRegions, clObj=clObj)
#'     # parallel::stopCluster(clObj)
#' }
#'
BaitCapture = function(sampleSheet, genome, baits, clObj=NULL){

  QuasRprj = GetQuasRprj(sampleSheet, genome)

  InBaits=QuasR::qCount(QuasRprj, baits, clObj = clObj)

  seq_length = seqlengths(genome)
  tiles = tileGenome(seq_length, tilewidth = max(seq_length), cut.last.tile.in.chrom=TRUE)
  GW = QuasR::qCount(QuasRprj, tiles, clObj=clObj)

  capture_efficiency = c()
  for(n in seq_along(unique(QuasRprj@alignments$SampleName))){
    capture_efficiency = c(capture_efficiency, sum(InBaits[,n+1]) / sum(GW[,n+1]))
  }

  return(capture_efficiency)

}

#' Subset Granges for given samples
#'
#' Inner utility for LowCoverageMethRateDistribution
#'
#' @param Single GRanges object as returned by CallContextMethylation function
#' @param Samples vector of sample names as they appear in the SampleName field of the QuasR sampleSheet
#'
#' @import GenomicRanges
#' @importFrom dplyr as_tibble
#'
SubsetGRangesForSamples = function(GRanges_obj, Samples){

  Cols2Discard = grep(paste0(c(Samples, "^GenomicContext$"), collapse = "|"), colnames(elementMetadata(GRanges_obj)), invert = TRUE)
  elementMetadata(GRanges_obj)[,Cols2Discard] = NULL
  All_NA_rows = is.na(rowMeans(as_tibble(elementMetadata(GRanges_obj)[,-ncol(elementMetadata(GRanges_obj))]), na.rm = TRUE))
  GRanges_obj_subset = GRanges_obj[!All_NA_rows]

  return(GRanges_obj_subset)

}

#' Manipulate GRanges into data.frame
#'
#' Inner utility for LowCoverageMethRateDistribution
#'
#' @param Single GRanges object as returned by CallContextMethylation function
#'
#' @import dplyr
#' @importFrom tidyr spread gather extract
#' @importFrom stats na.omit
#'
GRanges_to_DF = function(GRanges_obj){

  GRanges_obj %>%
    as_tibble() %>%
    select(matches("^seqnames$|^start$|_MethRate$|_Coverage$|^GenomicContext$|^Bins$")) %>%
    gather(Sample, Score, -seqnames, -start, -GenomicContext, -Bins) %>%
    na.omit() %>%
    extract(Sample, c("Sample", "Measurement"), regex = "(.*)_(Coverage$|MethRate|)") %>%
    spread(Measurement, Score) %>%
    mutate(Methylated = Coverage*MethRate) %>%
    group_by(Sample, GenomicContext, Bins) %>%
    summarise(TotCoverage = sum(Coverage), TotMethylated = sum(Methylated), ObservedMeth = TotMethylated/TotCoverage) %>%
    ungroup() %>%
    group_by(Sample, GenomicContext) %>%
    mutate(ExpectedMeth = seq(1,100,100/n())/100) %>%
    ungroup() %>%
    select(-TotCoverage, -TotMethylated, -Bins) -> DF

  return(DF)

}

#' Plot low coverage methylation rate
#'
#' Inner utility for LowCoverageMethRateDistribution
#'
#' @param Plotting_DF data.frame as returned by GRanges_to_DF function.
#'
#' @import ggplot2
#'
Plot_LowCoverageMethRate = function(Plotting_DF){

  ggplot(Plotting_DF, aes(x=ExpectedMeth,y=ObservedMeth, group=interaction(Sample,Coverage), color=Sample)) +
    geom_line(aes(linetype=Coverage, size=Coverage)) +
    ylab("Observed Methylation Rate") +
    xlab("Expected Methylation Rate") +
    facet_wrap(~GenomicContext) +
    theme_classic() +
    scale_size_manual("type", values = c(2.5, 1), guide = "none")

}

#' Low Coverage Methylation Rate MSE
#'
#' Calculate Mean squared error (MSE) of methylation rate distribution estimates for low coverage samples
#'
#' @param BinnedMethRate data.frame as returned by GRanges_to_DF function.
#'
#' @importFrom dplyr filter
#'
LowCoverageMethRate_MSE = function(BinnedMethRate){

  AllSamples = unique(BinnedMethRate$Sample)
  MSE_DF = data.frame(Sample = c(), MSE = c())
  for (sample in seq_along(AllSamples)){

    x = filter(BinnedMethRate, Sample == AllSamples[sample])$ExpectedMeth
    y = filter(BinnedMethRate, Sample == AllSamples[sample])$ObservedMeth

    MSE_DF = rbind(MSE_DF, data.frame(Sample = AllSamples[sample], MSE = (mean((y - x) ** 2))))

  }

  return(MSE_DF)

}

#' Plot Low Coverage Methylation Rate MSE
#'
#' Produce barplot of MSE values calculated for methylation rate distribution estimates of low coverage samples
#'
#' @param MSE_DF data.frame as returned by the LowCoverageMethRate_MSE function
#'
#' @import ggplot2
#'
Plot_LowCoverageMethRate_MSE = function(MSE_DF){

  MSE_DF %>%
    ggplot(aes(Sample, MSE)) +
    geom_bar(stat = 'identity') +
    ylab("Mean squared error\n(MSE)") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) -> pl

  return(pl)

}

#' Low Coverage Methylation Rate
#'
#' Monitor methylation rate distribution in a low coverage samples as compared to a high coverage "reference" one.
#' It bins cytosines with similar methylation rates (as observed in the HighCoverage sample) into bins. A single
#' methylation rate value is computed for each bin
#'
#' @param LowCoverage Single GRanges object as returned by CallContextMethylation function run with Coverage parameter set to 1. The object can also contain cytosines from multiple contexts
#' @param LowCoverage_samples Samples to use from the LowCoverage object. Either a string or a vector (for multiple samples).
#' @param HighCoverage Single GRanges object as returned by CallContextMethylation function. The object can also contain cytosines from multiple contexts.
#' @param HighCoverage_samples Single sample to use from HighCoverage. String
#' @param bins The number of bins for which to calculate the "binned" methylation rate. Defaults to 50
#' @param returnDF Whether to return the data.frame used for plotting. Defaults to FALSE
#' @param returnPlot Whether to return the plot. Defaults to TRUE
#' @param MSE Whether to calculate Mean squared error (MSE) of methylation rate distribution estimates for low coverage samples. Defaults to TRUE
#' @param return_MSE_DF  Whether to return a data.frame of computed MSE values. Defaults to FALSE
#' @param return_MSE_plot Whether to return a barplot of computed values. Defaults to TRUE
#'
#' @import GenomicRanges
#' @importFrom dplyr mutate
#'
#' @export
#'
#' @examples
#' # I don't have enough example data for this
#' # LowCoverageMethRateDistribution(LowCoverage = LowCoverage$DGCHN,
#' #                                 LowCoverage_samples = LowCoverage_Samples,
#' #                                 HighCoverage = HighCoverage$DGCHN,
#' #                                 HighCoverage_samples = HighCoverage_samples[1],
#' #                                 returnDF = FALSE,
#' #                                 returnPlot = TRUE,
#' #                                 MSE = TRUE,
#' #                                 return_MSE_DF = FALSE,
#' #                                 return_MSE_plot = TRUE)
#'
LowCoverageMethRateDistribution = function(LowCoverage, LowCoverage_samples, HighCoverage, HighCoverage_samples,
                                           bins = 50, returnDF = FALSE, returnPlot = TRUE,
                                           MSE = TRUE, return_MSE_DF = FALSE, return_MSE_plot = TRUE){

  message("Subsetting GRanges for given samples")
  HighCoverage = SubsetGRangesForSamples(HighCoverage, HighCoverage_samples)
  LowCoverage = SubsetGRangesForSamples(LowCoverage, LowCoverage_samples)

  message("Identifying bins on high coverage sample")
  HighCoverage_MethRate = elementMetadata(HighCoverage)[,grep("_MethRate$", colnames(elementMetadata(HighCoverage)))]
  HighCoverage$Bins = cut(HighCoverage_MethRate, breaks = seq(0,1,1/bins), include.lowest = TRUE)

  message("Assigning bins to low coverage samples")
  Overlaps = findOverlaps(LowCoverage, HighCoverage)
  LowCoverage$Bins = factor(NA, levels = levels(HighCoverage$Bins))
  LowCoverage$Bins[queryHits(Overlaps)] = HighCoverage$Bins[subjectHits(Overlaps)]

  message("Computing binned methylation rates")
  GRanges_objects = list(HighCoverage, LowCoverage)
  Coverages = c("High", "Low")
  BinnedMethRate = Reduce(rbind, lapply(seq_along(GRanges_objects), function(nr){
    mutate(GRanges_to_DF(GRanges_objects[[nr]]), Coverage = Coverages[nr])
  }))

  ReturnList = list()

  if (returnDF){
    ReturnList$MethylationDistribution_DF = BinnedMethRate
  }

  if (returnPlot) {
    Plot = Plot_LowCoverageMethRate(BinnedMethRate)
    ReturnList$MethylationDistribution_plot = Plot
  }

  if (MSE){

    MSE_DF = LowCoverageMethRate_MSE(BinnedMethRate)

    if (return_MSE_DF){
      ReturnList$MSE_DF = MSE_DF
    }

    if (return_MSE_plot){
      MSE_plot = Plot_LowCoverageMethRate_MSE(MSE_DF)
      ReturnList$MSE_plot = MSE_plot
    }

  }

  return(ReturnList)

}


# NEED TO DO: following function now doesn't work with external data, and correlation for example data is NA

#' Intersample correlation
#'
#' pair plot of sample correlations
#' @param samples Avg methylation object. Can also be set to "Example" to produce plot using example data of the kind specified by the @param CellType
#' @param context one of "AllCs", "DGCHN", "NWCGW". The first should be chosen for TKO experiments. For experiments carried on WT cells, we recommend checking both the "DGCHN" and "NWCGW" contexts by running this function once per context.
#' @param CellType Cell type to compare your samples to. At the moment, this can be one of "ES", "NP", "TKO".
#' @param saveAs Full path to output plot file
#'
#' # export
#'
#' # examples
#'
SampleCorrelation = function(samples, context, CellType, saveAs=NULL){

  # Get methylation data from previous SMF experiments
  # AllC = readRDS(system.file("extdata", "AllCreduced.rds", package = "SingleMoleculeFootprinting", mustWork = TRUE))
  # metMat_ref = readRDS(system.file("extdata", "ReducedRefMat.rds", package = "SingleMoleculeFootprinting", mustWork = TRUE))
  AllC = readRDS("/g/krebs/barzaghi/Rscripts/R_package/AllCreduced.rds")
  metMat_ref = readRDS("/g/krebs/barzaghi/Rscripts/R_package/ReducedRefMat.rds")

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
    # for (i in seq_along(ncol(contextMet))){
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
  comparison_mat[,seq_along(ncol(contextMet))] = contextMet
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
