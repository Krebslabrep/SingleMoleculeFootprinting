#author: arnaud.krebs@embl.de
#created 02.07.2020
#script to call methylation average in GC and CG context
#splits the call for individual chromosomes to avoid memory overload (only needed for large genomes)

library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(QuasR)
source('/g/krebs/barzaghi/Rscripts/R_package/SingleMoleculeFootprinting/R/functions/context_methylation_functions.r')

#############################
# Arguments
#############################
Qinput='../QuasR_input.txt' #QuasR alignement file containing a sample
out_path=paste("/g/krebs/barzaghi/Rscripts/R_package/")
cO=19 #minimal coverage
nb.cores=20 #nb of cores to be used

# Create QuasR project
NOMEproj=qAlign(sampleFile=Qinput,
                genome="BSgenome.Mmusculus.UCSC.mm10",
                projectName = "CTCF_amplicon",
                paired="fr",
                bisulfite="undir")
NOMEaln=as.data.frame(alignments(NOMEproj)[[1]])
NOMEproj@aligner = "Rbowtie"

#############################
# 1. Single site example
#############################
Region_of_interest = GRanges(seqnames = "chr1", ranges = IRanges(start = 31210117, end = 31210616), strand = "+")

# Quantify methylation in region of interest
meth_gr <- qMeth(NOMEproj, mode="allC", Region_of_interest)
# Estract methylation info from contexts of interest only, the function also collapses strands and filters for coverage
contextMet=CallContextMethylation(meth_gr, cO, Mmusculus)

#############################
# 2. Genome-wide example:
# given the example data only cover one genomic region,
# most of the of the output files will contain no methylation information
#
# This process can be quite lengthy and computationally demanding, it is advisable to run it on the cluster
#############################

# Partition a genome by chromosome ("natural partitioning")
musmus_length=seqlengths(Mmusculus)
tiles <- tileGenome(musmus_length, tilewidth=max(musmus_length),cut.last.tile.in.chrom=TRUE)

# Call the methylation genome wide for all Cs, loop/chromosome
cluObj=makeCluster(nb.cores)
lapply(1:21,function(i){

  print(i)
	meth_gr <- qMeth(NOMEproj, mode="allC", tiles[1:21][i], clObj=cluObj)
	contextMet=CallContextMethylation(meth_gr, cO, Mmusculus)
	saveRDS(contextMet, paste0(out_path,'/Context_met_call_',NOMEproj@projectName,'_',as.character(seqnames( tiles[1:21][i])),'_Co',as.character(cO),'.rds',sep=''))

	})

# Filter away Cytosines with low coverage in all samples and combine chromosome-specific objects
AllCf=mclapply(1:21,function(i){

	contextMet=readRDS(paste(out_path,'/Context_met_call_',NOMEproj@projectName,'_',as.character(seqnames( tiles[1:21][i])),'_Co',as.character(cO),'.rds',sep=''))
	CG=contextMet[[1]]
	GC=contextMet[[2]]
	AllC=c(CG,GC)
	met=elementMetadata(AllC)
	met2=met[,1:(ncol(met)-1)]
	cov.inx=!rowSums(is.na(met2))==ncol(met2)
	AllCf=AllC[cov.inx]
	AllCf

}, mc.cores=nb.cores)

AllC=unlist(GRangesList(AllCf))
AllC=sort(AllC)

# save final object
saveRDS(AllC, paste0(out_path,'/Context_methylation_call_',NOMEproj@projectName,'.rds'))

# remove chromosome-wise temporary files
lapply(1:21, function(i){
  file.remove(paste0(out_path,'/Context_met_call_',NOMEproj@projectName,'_',as.character(seqnames( tiles[1:21][i])),'_Co',as.character(cO),'.rds',sep=''))
})

