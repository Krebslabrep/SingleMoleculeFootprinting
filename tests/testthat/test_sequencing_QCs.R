library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)

CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
Qinput = paste0(CacheDir, "/NRF1Pair_Qinput.txt")

test_that("BaitCapture returns a numeric", {

  expect_true(is.numeric(BaitCapture(sampleSheet = Qinput,
                                     genome = BSgenome.Mmusculus.UCSC.mm10,
                                     baits = EnrichmentRegions_mm10.rds()))
    )

})
