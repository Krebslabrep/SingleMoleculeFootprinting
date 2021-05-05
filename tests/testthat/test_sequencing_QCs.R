library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)

test_that("BaitCapture returns a numeric", {

  expect_true(is.numeric(BaitCapture(sampleSheet = paste0(CacheDir, "/NRF1Pair_Qinput.txt"),
                                     genome = BSgenome.Mmusculus.UCSC.mm10,
                                     baits = EnrichmentRegions_mm10.rds()))
    )

})
