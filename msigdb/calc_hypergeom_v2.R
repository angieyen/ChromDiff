# Copyright 2015 Angela Yen
# This file is part of ChromDiff.
# ChromDiff is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ChromDiff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ChromDiff.  If not, see <http://www.gnu.org/licenses/>.
## reldir should be relative path to main ChromDiff_v2 directory. Should not need to be changed unless msigdb directory is moved
reldir="../"
source(paste(reldir, "setvars.R", sep=""), chdir=T)
testone=FALSE
if(testone) {
	property="specialgi"
	a.option="BRAIN"
	b.option="GI"
	test="wilcox"
	correction="fdr"
	model="core_gencode_v10"
	metadatafile=paste0(reldir, "data/final_celltype_metadata.txt")
	genefile=paste0(reldir, "data/gencode_genes_full.txt")
	covariate_mat_file=paste0(reldir, "data/cov.mat.txt")
	map_covariates_file=paste0(reldir, "data/map_vars_covariates.txt")
	expfile=paste0(reldir, "rnaseq/57epigenomes.RPKM.pc")
	state_annotations_file=paste0(reldir, "data/core_annotation.txt")
	generegions_label="gencode_v10"
}else{
	args<-commandArgs(TRUE)
	property=args[1]
	a.option=args[2]
	b.option=args[3]
	test=args[4]
	correction=args[5]
	model=args[6]
    metadatafile=args[7]
    genefile=args[8]
    covariate_mat_file=args[9]
    map_covariates_file=args[10]
    expfile=args[11]
    state_annotations_file=args[12]
	generegions_label=args[13]
}
a.label=a.option
b.label=b.option
metric="perc"
clusts=TRUE
overwrite=TRUE
idtype="symbols"
msigdbset="c2_c5"
plottype="domstate"
suffix="_matched"
set_variables(model, metadatafile,  genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label)
source(paste(reldir, "funcs.R", sep=""), chdir=T)
source(paste(reldir, "clust_funcs.R", sep=""), chdir=T)
source("msigdb_funcs.R")

## load MSigDB gene sets
msigdata=get_msigdata(msigdbset, idtype)
plottypesuffix=get_plottype_suffix(plottype)
plotdir=paste0(reldir,get_full_plotdir(model, a.label, b.label, test, metric, correction))
## read in gencode ids, convert to ensembl ids then entrez ids

allpvalfile=paste0(reldir, get_pval_file(model, metric, test, correction, property, a.label, b.label))
if(file.exists(allpvalfile)) {
	outfile=get_msigdb_resultsfile(reldir, msigdbset, idtype, model,  metric, test, correction, property, a.label, b.label)
	dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
	if((!overwrite) && file.exists(outfile) && file.info(outfile)$size>0) {
		quit()	
	}
	print(paste("Calculating enrichments for ", plotdir, "...", sep=""))
	pval.table=read.table(allpvalfile)
	names(pval.table)=c("feat", "pval")
	if(nrow(pval.table)!=length(get_featnames(model, metric, reldir))) {
		stop(paste("Problem calculating enrichments for ", plotdir, ": wrong number of entries", sep=""))
	}
	gencode.genes=get_sig_genes(pval.table)
	if(length(gencode.genes)>0) {
		matched.genes=convert(gencode.genes, idtype)
		final.genes=remove.extra(matched.genes)

		print("Getting enrichments...")
		## get overall msig enrichments
		gene.universe=get_gene.universe(idtype)
		results=get_enrichments(msigdata, final.genes, gene.universe)

		## if matrix  and has more than 0 rows, or if it is a vector (has 1 result)
		if(!is.matrix(results) || nrow(results)>0) {
			dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
			write.table(results, quote=FALSE, file=outfile, row.names=FALSE, sep="\t")
		}
		rm(results)

		if(clusts) {
			vert.ints.file=get_vert.ints.file(plotdir, plottypesuffix)
			if(!file.exists(vert.ints.file)) {
				stop(paste0("File ", vert.ints.file, " does not exist.")) 
			}
			load(vert.ints.file)
			if(!(exists("vert.ints.noreps") && exists("minclustsize"))) {
				stop(paste0(vert.ints.file, " did not load vert.ints.noreps and minclustsize."))
			}
			if(!(is.list(vert.ints.noreps) && length(vert.ints.noreps)>0)) {
				stop(paste0("vert.ints.noreps was not a list or was length 0"))
			}
			vert.ints.noreps=order.list.by.first(vert.ints.noreps)
			## load gene order to match plot and clusterings
			geneorderfile=get_geneorder_file(plotdir, suffix, plottypesuffix)
			if(!file.exists(geneorderfile)) {
				stop(paste0(geneorderfile, " does not exist, so can not run clustering analysis."))
			}
			gencode.genes=as.vector(t(read.table(geneorderfile, skip=1)))
			final.genes=convert(gencode.genes, idtype)
			clustcount=1
			for(clustind in 1:length(vert.ints.noreps)) {
				currclust.indices=vert.ints.noreps[[clustind]]
				if((currclust.indices[2]-currclust.indices[1])>minclustsize) {
					outfile=get_msigdb_clustfile(reldir, msigdbset, idtype, model, metric, test, correction, property, a.label, b.label, clustcount)
					if(!overwrite && file.exists(outfile) && file.info(outfile)$size>0) {
						next
					}
					## get genes in current cluster
					currclust=unique(final.genes[(currclust.indices[1]-0.5):(currclust.indices[2]+.5)])
					finalclust=remove.extra(currclust)
					results=get_enrichments(msigdata, finalclust, gene.universe)
					if(!is.matrix(results) || nrow(results)>0) {
						dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
						write.table(results, quote=FALSE, file=outfile, row.names=FALSE, sep="\t")
					}
					rm(results)
					clustcount=clustcount+1
				}
			}
		}
	}
}
