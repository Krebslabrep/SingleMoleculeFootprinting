library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(testthat)

CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
Qinput = paste0(CacheDir, "/NRF1Pair_Qinput.txt")
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
Region_of_interest <- GRanges(seqnames = "chr6",
                              ranges = IRanges(start = 88106000, end = 88106500),
                              strand = "*")
Methylation <- CallContextMethylation(sampleSheet = Qinput,
                                      sample = MySample,
                                      genome = BSgenome.Mmusculus.UCSC.mm10,
                                      range = Region_of_interest,
                                      coverage = 20,
                                      ConvRate.thr = 0.2)
Meth_SM = Methylation[[2]]
TFBS <- GenomicRanges::GRanges("chr6",
                                IRanges(88106253, 88106263),
                                strand = "-")

test_that("BinMethylation returns named numeric vector", {

  expect_true(is.numeric(BinMethylation(MethSM = Meth_SM, TFBS = TFBS, bin = c(-15,15))))
  expect_true(!is.null(names(BinMethylation(MethSM = Meth_SM, TFBS = TFBS, bin = c(-15,15)))))

})

test_that("SortReads returns a list", {

  expect_true(is.list(SortReads(MethSM = Meth_SM, TFBS = TFBS, list(c(-35,-25), c(-15,15), c(25,35)), SortByCluster = FALSE)))

})
