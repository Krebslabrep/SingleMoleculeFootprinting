<div style="text-align: justify;">

# SingleMoleculeFootprinting

# Introduction
*SingleMoleculeFootprinting* is an R package build around [QuasR](https://github.com/fmicompbio/QuasR) and tailored for the analysis Single Molecule Footprinting (SMF) data.

SMF is a high-throughput sequencing technology developed in the [Krebs laboratory](https://www.embl.de/research/units/genome_biology/krebs/index.html). It consists of marking accessible genomic cytosines in the GpC and CpG contexts using exogenous methytransferase enzymes and of subsequently performing bisulfite sequencing (BS). Consequently, cytosines protected by the binding of DNA-interacting proteins (e.g. TFs, nucleosomes, GTFs, etc..) will result unmethylated, while the accessible cytosines will result methylated.

With the present package, we provide functions to perform basic SMF data analysis starting from aligned bam files up to the biological interpretation of results over single sites.

# Preamble: preprocessing SMF data
To ensure compatibility with our downstream tools, we recommend aligning sequencing reads using the QuasR function [qAlign](https://www.rdocumentation.org/packages/QuasR/versions/1.12.0/topics/qAlign) as follows
```r
cl = makeCluster(40)
prj = QuasR::qAlign(sampleFile = sampleFile,
              genome = genome,
              aligner = "Rbowtie",
              projectName = "prj", 
              paired = "fr",
              bisulfite = "undir", 
              alignmentParameter = "-e 70 -X 1000 -k 2 --best -strata",
              alignmentsDir = "./", 
              cacheDir = tempdir(),
              clObj = cl)
```
For more details on how to structure the **sampleFile** argument we refer to the [qAlign](https://www.rdocumentation.org/packages/QuasR/versions/1.12.0/topics/qAlign) documentation.
For more details on SMF data preprocessing we refer to the computational steps of our methods manuscript *link to pre-print*.

# Installation
To install *SingleMoleculeFootprinting*, execute the following
```r
remotes::install_github(repo = "https://github.com/Krebslabrep/SingleMoleculeFootprinting.git", ref = "main", build_vignettes = FALSE)
```

# SingleMoleculeFootprinting usage
For instructions on how to use the *SingleMoleculeFootprinting* package and an example of analysis, consult our [vignette](https://htmlpreview.github.io/?https://github.com/Krebslabrep/SingleMoleculeFootprinting/blob/main/vignettes/SingleMoleculeFootprinting.html)


</div>
