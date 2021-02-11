#' Get QuasRprj
#'
#' @param sampleSheet QuasR pointer file
#' @param genome BSgenome
#'
#' @import QuasR
#'
#' @export
#'
GetQuasRprj = function(sampleSheet, genome){

  QuasRprj=QuasR::qAlign(sampleFile=sampleSheet,
                        genome=genome@pkgname,
                        projectName = "NRF1pair_DE_example",
                        paired="fr",
                        aligner = "Rbowtie",
                        bisulfite="undir")
  QuasRprj@aligner = "Rbowtie"

  return(QuasRprj)

}

#' Get Single Molecule methylation matrix
#'
#' Used internally as the first step in CallContextMethylation
#'
#' @param QuasRprj QuasR project object as returned by calling the QuasR function qAling on previously aligned data
#' @param range GenimocRange representing the genomic region of interest
#' @param sample One of the sample names as reported in the SampleName field of the QuasR pointer file provided to qAlign. N.b. all the files with the passed sample name will be used to call methylation
#'
#' @import QuasR
#' @importFrom data.table data.table
#' @importFrom data.table dcast
#'
#' @return Single molecule methylation matrix (all Cytosines)
#'
#' @export
#'
GetSingleMolMethMat<-function(QuasRprj,range,sample){

  Cs=QuasR::qMeth(QuasRprj[grep(sample, QuasRprj@alignments$SampleName)], query=range, mode="allC",reportLevel="alignment", collapseBySample = TRUE)
  # use data.table to get a 1,0 matrix of methylation profiles
  # all.cids=unique(Cs[[sample]]$Cid) # get all possible C locations
  # make the data.table object
  dt=data.table::data.table(meth=Cs[[sample]]$meth, aid=Cs[[sample]]$aid, cid=Cs[[sample]]$Cid)
  spread_dt = data.table::dcast(dt, aid ~ cid, value.var = "meth") ## MEMORY issue with GW data
  meth_matrix = as.matrix(spread_dt, rownames = "aid")
  return(meth_matrix)

}

#' Calculate reads conversion rate
#'
#' @param MethSM as comes out of the func GetSingleMolMethMat
#' @param chr Chromosome, MethSM doesn't carry this info
#' @param genome BSgenome
#' @param thr Double between 0 and 1. Threshold above which to filter reads. Defaults to 0.2
#'
#' @import GenomicRanges
#' @import Biostrings
#' @importFrom IRanges IRanges
#'
#' @return Filtered MethSM
#'
FilterByConversionRate = function(MethSM, chr, genome, thr=0.2){

  CytosineRanges = GRanges(chr,IRanges(as.numeric(colnames(MethSM)),width = 1))
  GenomicContext = Biostrings::getSeq(genome, resize(CytosineRanges,3,fix='center'))
  IsInContext = Biostrings::vcountPattern('GC',GenomicContext)==1 | vcountPattern('CG',GenomicContext)==1

  # filter based on conversion
  ConvRate=rowMeans(MethSM[,!IsInContext],na.rm=TRUE)
  message(paste0(round((sum(ConvRate>=thr)/length(ConvRate))*100, digits = 2), "% of reads found with conversion rate above ", thr))
  FilteredSM = MethSM[ConvRate<thr,]

  return(FilteredSM)

}

#' Detect type of experiment
#'
#' @param Samples SampleNames field from QuasR sampleSheet
#'
DetectExperimentType = function(Samples){

  if(length(grep("_NO_", Samples)) > 0){
    ExpType = "NO"
  }else if(length(grep("_SS_",Samples)) > 0){
    ExpType = "SS"
  }else if(length(grep("_DE_",Samples)) > 0){
    ExpType = "DE"
  }else{stop("Sample name error. No _McvPI_ , _SsssI_ or _McvPISsssI_")}

  message(paste0("Detected experiment type: ", ExpType))
  return(ExpType)

}

#' Filter Cytosines in context
#'
#' @param MethGR Granges obj of average methylation
#' @param genome BSgenome
#' @param context Context of interest (e.g. "GC", "CG",..)
#'
#' @import GenomicRanges
#' @import Biostrings
#'
#' @return filtered Granges obj
#'
FilterContextCytosines <- function(MethGR, genome, context){

  CytosineRanges = GRanges(seqnames(MethGR), ranges(MethGR), strand(MethGR)) # performing operation without metadata to make it lighter
  seqinfo(CytosineRanges) <- seqinfo(MethGR)
  GenomicContext = Biostrings::getSeq(genome, suppressWarnings(trim(IRanges::resize(CytosineRanges, width = 5, fix = "center"))), as.character = FALSE) # I checked: it is strand-aware
  
  # Fix trimming issues necessary for Genome-Wide analysis
  trimmed <- which(!width(GenomicContext) == 5)
  for(k in trimmed){
    if(as.character(strand(CytosineRanges[k])) == '+'){
      GenomicContext[k] <- paste0(c(rep('N', 5-width(GenomicContext[k])), as.character(GenomicContext[k])), collapse="")
    } else if(as.character(strand(CytosineRanges[k])) == '-'){
      GenomicContext[k] <- paste0(c(as.character(GenomicContext[k]), rep('N', 5-width(GenomicContext[k]))), collapse="")
    }
  }

  IsInContext = unlist(lapply(Biostrings::vmatchPattern(context, subject = GenomicContext, fixed = FALSE), function(i){length(i)>0}))

  # elementMetadata(CytosineRanges)$GenomicContext = GenomicContext
  # elementMetadata(CytosineRanges)$IsInContext = IsInContext
  MethGR_InContext = MethGR[IsInContext]

  return(MethGR_InContext)

}

#' Fixing overhang before stand collapsing
#'
#' @param MethGR Granges obj of average methylation
#' @param context context
#' @param which "Top"|"Bottom"
#'
#' @import GenomicRanges
#' @import BiocGenerics
#' @importFrom IRanges IRanges
#'
#' @return MethGR with fixed overhang
#'
FixOverhang = function(MethGR, context, which){

  if (which == "Top"){
    newParams = list(min(start(MethGR)) - 1, ifelse(context == "CG", "+", "-"))
  } else if (which == "Bottom"){
    newParams = list(max(start(MethGR)) + 1, ifelse(context == "CG", "-", "+"))
  }

  overhang_fix = GRanges(seqnames = unique(as.character(seqnames(MethGR))), ranges = IRanges(newParams[[1]], width = 1), strand = newParams[[2]])
  fix_metadata = data.frame(matrix(0, ncol = length(colnames(elementMetadata(MethGR))), nrow = 1))
  colnames(fix_metadata) = colnames(elementMetadata(MethGR))
  elementMetadata(overhang_fix) = fix_metadata

  if (which == "Top"){
    MethGR = append(overhang_fix, MethGR)
  } else if (which == "Bottom"){
    MethGR = append(MethGR, overhang_fix)
  }

  return(MethGR)

}

#' Collapse strands
#'
#' @param MethGR Granges obj of average methylation
#' @param context "GC" or "CG". Broad because indicates just the directionality of collapse.
#'
#' @import GenomicRanges
#'
#' @return MethGR with collapsed strands (everything turned to - strand)
#'
CollapseStrands = function(MethGR, context){

  TopStrandToFix = ifelse(context == "CG", "-", "+") # if this is GR first strand, we need to fix
  # GC needs to start with "-", CG with "+"
  if (as.character(strand(MethGR[1]))==TopStrandToFix){
    message("Strand collapsing: Fixing top (left) overhang")
    MethGR = suppressWarnings(FixOverhang(MethGR, context, "Top"))
  }
  # CG needs to start with "+", CG with "-"
  BottomStrandToFix = ifelse(context == "CG", "+", "-") # if this is GR last strand, we need to fix
  if (as.character(strand(MethGR[length(MethGR)]))==BottomStrandToFix){
    message("Strand collapsing: Fixing bottom (right) overhang")
    MethGR = suppressWarnings(FixOverhang(MethGR, context, "Bottom"))
  }

  # Initiating collapsed GR
  MethGR_collapsed <- MethGR[seq(ifelse(context == "CG", 2, 1),length(MethGR),by=2)]
  # Resizing
  # if (context == "CG"){
  #   end(MethGR_collapsed) <- end(MethGR_collapsed)+1
  # } else if (context == "GC"){
  #   start(MethGR_collapsed) <- start(MethGR_collapsed)-1
  # }
  # Summing the + and - coverages
  values(MethGR_collapsed) <- as.matrix(values(MethGR[seq(1,length(MethGR),by=2)])) + as.matrix(values(MethGR[seq(2,length(MethGR),by=2)]))

  return(MethGR_collapsed)

}

#' Collapse strands in single molecule matrix
#'
#' The idea here is that (regardless of context) if a C is on the - strand, calling getSeq on that coord (N.b. unstranded, that's the important bit) will give a "G', a "C" if it's a + strand.
#'
#' @param MethSM Single molecule matrix
#' @param context "GC" or "CG". Broad because indicates just the directionality of collapse.
#' @param genome BSgenome
#' @param chr Chromosome, MethSM doesn't carry this info
#'
#' @import GenomicRanges
#' @import Biostrings
#' @importFrom plyr rbind.fill.matrix
#' @importFrom IRanges IRanges
#'
#' @return Strand collapsed MethSM
CollapseStrandsSM = function(MethSM, context, genome, chr){

  CytosineRanges = GRanges(chr,IRanges(as.numeric(colnames(MethSM)),width = 1))
  GenomicContext = Biostrings::getSeq(genome, CytosineRanges)
  IsMinusStrand = GenomicContext=="G" # no need to select for GC/CG context --> MethSM passed already subset

  # Separate MethSM into - and + strands
  MethSM_minus = MethSM[apply(MethSM[,IsMinusStrand], 1, function(i){sum(is.na(i)) != length(i)}) > 0, IsMinusStrand]
  MethSM_plus = MethSM[apply(MethSM[,!IsMinusStrand], 1, function(i){sum(is.na(i)) != length(i)}) > 0, !IsMinusStrand]
  NrPlusReads = dim(MethSM_plus[1])
  message(paste0(ifelse(is.null(NrPlusReads), 0, NrPlusReads), " reads found mapping to the + strand, collapsing to -"))

  # Turn + into -
  offset = ifelse(context == "GC", -1, +1) # the opposite if I was to turn - into +
  colnames(MethSM_plus) = as.character(as.numeric(colnames(MethSM_plus)) + offset)

  # Merge matrixes
  StrandCollapsedSM = plyr::rbind.fill.matrix(MethSM_minus, MethSM_plus)
  rownames(StrandCollapsedSM) = c(rownames(MethSM_minus), rownames(MethSM_plus))
  StrandCollapsedSM = StrandCollapsedSM[,as.character(sort(as.numeric(colnames(StrandCollapsedSM))))]

  return(StrandCollapsedSM)

}

#' Filter Cs for coverage
#'
#' @param MethGR Granges obj of average methylation
#' @param thr converage threshold
#'
#' @import GenomicRanges
#'
#' @return filtered MethGR
#'
CoverageFilter <- function(MethGR, thr){

  Tcounts = elementMetadata(MethGR)[seq(1, ncol(elementMetadata(MethGR)), by = 2)]
  Mcounts = elementMetadata(MethGR)[seq(2, ncol(elementMetadata(MethGR)), by = 2)]
  MethylationMatrix = as.matrix(Mcounts)/as.matrix(Tcounts)

  #filter for coverage
  CovFilter=as.matrix(Tcounts)>thr
  for (i in seq_along(ncol(MethylationMatrix))){
    MethylationMatrix[!CovFilter[,i],i] = NA
  }

  #bind the GRanges with the scores
  elementMetadata(MethGR) = NULL
  elementMetadata(MethGR) = MethylationMatrix

  # do the actual filtering
  MethGRFiltered = MethGR[!(rowSums(is.na(as.matrix(elementMetadata(MethGR)))) == length(elementMetadata(MethGR)))]

  return(MethGRFiltered)

}

#' Call Context Methylation
#'
#' Can deal with multiple samples
#'
#' @export
#'
#' @param sampleSheet QuasR pointer file
#' @param sample for now this works for sure on one sample at the time only
#' @param genome BSgenome
#' @param range GenimocRange representing the genomic region of interest
#' @param coverage coverage threshold. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to 0.2
#'
#' @import QuasR
#' @import GenomicRanges
#' @import BiocGenerics
#'
#' @return List with two Granges objects: GC and CG average methylation
#'
CallContextMethylation=function(sampleSheet, sample, genome, range, coverage = 20, ConvRate.thr = 0.2, singleMolecules = TRUE){

  message("Setting QuasR project")
  QuasRprj = GetQuasRprj(sampleSheet, genome)
  Samples = QuasR::alignments(QuasRprj)[[1]]$SampleName

  message("Calling methylation at all Cytosines")
  MethGR = QuasR::qMeth(QuasRprj[grep(sample, Samples)], mode="allC", range, collapseBySample = TRUE, keepZero = TRUE)
  if(singleMolecules){
    MethSM = GetSingleMolMethMat(QuasRprj, range, sample) # this selects the sample internally ---> TO FIX
    MethSM = FilterByConversionRate(MethSM, chr = seqnames(range), genome = genome, thr = ConvRate.thr)
  }

  message("Subsetting Cytosines by permissive genomic context (NGCNN, NNCGN)") # Here we use a permissive context: needed for the strand collapsing
  ContextFilteredMethGR = list(GC = FilterContextCytosines(MethGR, genome, "NGCNN"),
                               CG = FilterContextCytosines(MethGR, genome, "NNCGN"))
  if(singleMolecules){
    ContextFilteredMethSM = lapply(seq_along(ContextFilteredMethGR), function(i){MethSM[,colnames(MethSM) %in% as.character(start(ContextFilteredMethGR[[i]]))]})
  }

  message("Collapsing strands")
  StrandCollapsedMethGR = list(GC = CollapseStrands(ContextFilteredMethGR[[1]], context = "GC"),
                               CG = CollapseStrands(ContextFilteredMethGR[[2]], context = "CG"))
  if(singleMolecules){
    StrandCollapsedMethSM = list(GC = CollapseStrandsSM(ContextFilteredMethSM[[1]], context = "GC", genome = genome, chr = as.character(seqnames(range))),
                                 CG = CollapseStrandsSM(ContextFilteredMethSM[[2]], context = "CG", genome = genome, chr = as.character(seqnames(range))))
  }

  message("Filtering Cs for coverage")
  CoverageFilteredMethGR = list(GC = CoverageFilter(StrandCollapsedMethGR[[1]], thr = coverage),
                                CG = CoverageFilter(StrandCollapsedMethGR[[2]], thr = coverage))
  if(singleMolecules){
    CoverageFilteredMethSM = lapply(seq_along(CoverageFilteredMethGR), function(i){StrandCollapsedMethSM[[i]][,colnames(StrandCollapsedMethSM[[i]]) %in% as.character(start(CoverageFilteredMethGR[[i]]))]})
  }

  # Determining stric context based on ExpType
  ExpType = DetectExperimentType(Samples)
  if (ExpType == "NO"){
    ExpType_contexts = c("DGCHN", "NWCGW")
  } else if (ExpType == "SS"){
    ExpType_contexts = c("", "CG")
  } else if (ExpType == "DE"){
    ExpType_contexts = c("GCH", "HCG")
  }

  message(paste0("Subsetting Cytosines by strict genomic context (", ExpType_contexts[1], ", ", ExpType_contexts[2],") based on the detected experiment type: ", ExpType))
  ContextFilteredMethGR_strict = list(FilterContextCytosines(CoverageFilteredMethGR[[1]], genome, ExpType_contexts[1]),
                                      FilterContextCytosines(CoverageFilteredMethGR[[2]], genome, ExpType_contexts[2]))
  names(ContextFilteredMethGR_strict) = ExpType_contexts
  if(singleMolecules){
    ContextFilteredMethSM_strict = lapply(seq_along(ContextFilteredMethGR_strict), function(i){CoverageFilteredMethSM[[i]][,colnames(CoverageFilteredMethSM[[i]]) %in% as.character(start(ContextFilteredMethGR_strict[[i]]))]})
  }

  if (ExpType == "DE"){
    message("Merging matrixes")
    MergedGR = sort(append(ContextFilteredMethGR_strict[[1]], ContextFilteredMethGR_strict[[2]]))

    if(singleMolecules){
      # When some reads only cover either DGCHN or NWCGW positions cbind complains
      if(nrow(ContextFilteredMethSM_strict[[1]]) != nrow(ContextFilteredMethSM_strict[[2]]) | !(all(sort(rownames(ContextFilteredMethSM_strict[[1]])) %in% sort(rownames(ContextFilteredMethSM_strict[[2]]))))){

        # Which reads cover only one context?
        DGCHNonly_reads = rownames((ContextFilteredMethSM_strict[[1]]))[!(rownames((ContextFilteredMethSM_strict[[1]])) %in% rownames((ContextFilteredMethSM_strict[[2]])))]
        NWCGWonly_reads = rownames((ContextFilteredMethSM_strict[[2]]))[!(rownames((ContextFilteredMethSM_strict[[2]])) %in% rownames((ContextFilteredMethSM_strict[[1]])))]
        # Fill two dummy matrices to make the reads equal in the ContextFilteredMethSM_strict matrices
        DGCHNonly_mat = matrix(data = NA, nrow = length(DGCHNonly_reads), ncol = ncol(ContextFilteredMethSM_strict[[2]]), dimnames = list(DGCHNonly_reads, colnames(ContextFilteredMethSM_strict[[2]])))
        NWCGWonly_mat = matrix(data = NA, nrow = length(NWCGWonly_reads), ncol = ncol(ContextFilteredMethSM_strict[[1]]), dimnames = list(NWCGWonly_reads, colnames(ContextFilteredMethSM_strict[[1]])))
        # merge the ContextFilteredMethSM_strict matrices to the respective dummy
        ContextFilteredMethSM_strict[[1]] = BiocGenerics::rbind(ContextFilteredMethSM_strict[[1]], NWCGWonly_mat)
        ContextFilteredMethSM_strict[[2]] = BiocGenerics::rbind(ContextFilteredMethSM_strict[[2]], DGCHNonly_mat)

      }

      # Sort reads alphanumerically before binding, because cbind doesn't join (I've tested) matrices by rownames
      ContextFilteredMethSM_strict[[1]] = ContextFilteredMethSM_strict[[1]][sort(rownames(ContextFilteredMethSM_strict[[1]])),]
      ContextFilteredMethSM_strict[[2]] = ContextFilteredMethSM_strict[[2]][sort(rownames(ContextFilteredMethSM_strict[[2]])),]
      MergedSM = BiocGenerics::cbind(ContextFilteredMethSM_strict[[1]], ContextFilteredMethSM_strict[[2]])
      MergedSM = MergedSM[,as.character(sort(as.numeric(colnames(MergedSM))))]
    } else {
      MergedSM = NULL
    }
  }

  return(list(MergedGR, MergedSM))

}
