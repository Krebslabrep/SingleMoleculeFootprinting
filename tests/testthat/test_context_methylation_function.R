library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(testthat)

CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
Qinput = paste0(CacheDir, "/NRF1Pair_Qinput.txt")
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
Region_of_interest <- GRanges(seqnames = "chr6",
                              ranges = IRanges(start = 88106000, end = 88106500),
                              strand = "*")
QuasRprj = GetQuasRprj(Qinput, BSgenome.Mmusculus.UCSC.mm10)

test_that("GetSingleMolMethMat returns a matrix", {

  expect_true(is.matrix(GetSingleMolMethMat(QuasRprj, range = Region_of_interest, sample = MySample)))

})

test_that("DetectExperimentType behaves as expected", {

  Samples = QuasR::alignments(QuasRprj)[[1]]$SampleName
  expect_true(is.character(DetectExperimentType(Samples = Samples)))
  expect_message(DetectExperimentType(Samples = Samples, verbose = TRUE))

})


