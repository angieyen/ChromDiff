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

source("setvars.R", chdir=T)
source("mkplots.R", chdir=T)
source("funcs.R", chdir=T)
source("heatmap.3.R", chdir=T)
test=FALSE
if(test) {
	metric="perc"
	testtype="wilcox"
	correction="fdr"
	property="sex"
	group1="Female"
	group2="Male"
	label1="Female"
	label2="Male"
	model="core_gencode_v10"
	plottype="domstate"
	metadatafile="data/final_celltype_metadata.txt"
	genefile="data/gencode_genes_1kb_up.txt"
	covariate_mat_file="data/cov.mat.txt"
	map_covariates_file="data/map_vars_covariates.txt"
	expfile="data/57epigenomes.RPKM.pc"
	state_annotations_file="data/core_annotation.txt"
	generegions_label="gencode_v10_1kb"
} else {
	args<-commandArgs(TRUE)
	metric=args[1]
	testtype=args[2]
	correction=args[3]
	property=args[4]
	group1=args[5]
	group2=args[6]
	label1=args[7]
	label2=args[8]
	model=args[9]
	plottype=args[10]
    metadatafile=args[11]
    genefile=args[12]
    covariate_mat_file=args[13]
    map_covariates_file=args[14]
    expfile=args[15]
    state_annotations_file=args[16]
	generegions_label=args[17]
}
set_variables(model, metadatafile,  genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label)

convert2ensembl <- function(gencode.ids) {
	split.mat=sapply(gencode.ids, function(x) {return(strsplit(x, "[.]")[[1]])})
	return(split.mat[1,])
}
useIds=FALSE
useNames=FALSE
suffix2=""
suffix="_matched"	
plottypesuffix=get_plottype_suffix(plottype)
list[reorder_by_assocstate, useDend, reorder_by_domstate] = get_plottype_bools(plottype)

dend=FALSE
dend.choice="none"

plotdir=get_full_plotdir(model, label1, label2, testtype, metric, correction)
plotfile=paste(plotdir, "sig_exp_plot", suffix, suffix2, plottypesuffix,".pdf", sep="")

print("Loading in celltype, gene, and chromatin state ordering...")
## load in celltype ordering
celltypeorder.file=paste(plotdir, "sig_maj_celltypeorder", suffix, plottypesuffix, ".txt", sep="")
celltype.order=as.vector(t(read.table(celltypeorder.file, header=TRUE)))

## load in saved gene order and associated chromatin states
geneorder.file=paste(plotdir, "sig_maj_geneorder", suffix, plottypesuffix, ".txt", sep="")
geneorder=as.vector(t(read.table(geneorder.file, header=TRUE)))
## CONVERT original gene ids to ensembl gene ids
ensembl.geneorder=convert2ensembl(geneorder)	

## if no features, return, if only one feature, special plotting
if(length(ensembl.geneorder)==0) {
	print(paste( "No genes for", plotfile))
	quit()
}

chrstateorder.file=paste(plotdir, "sig_maj_chrstateorder", suffix, plottypesuffix, ".txt", sep="")
mapgenes2chrstates=as.vector(t(read.table(chrstateorder.file, header=TRUE)))

mat.to.plot=get_logexp(celltype.order, ensembl.geneorder)
if(reorder_by_domstate) {
	majdendfile=paste(plotdir, "sig_maj_plot_dend", plottypesuffix, ".rdata", sep="")
	#hclust=hclust(dist(t(mat.to.plot), method="euclidean"), method="complete")
	dend=load_or_get_dend(majdendfile, mat.to.plot, 2, TRUE,"dend")
	#dend=as.dendrogram(hclust)
	dend.choice="column"
}
if(!is.null(mat.to.plot)) {
	print("Generating final heatmap...")
	library(gplots)

	## pull out genes and celltypes we want 
	## make row colors that show relevant grouping
	pair.label=get_pair_label(label1, label2)
    celltypes.list=get_valid_celltypes(metric, property, pair.label)
	list1=celltypes.list$a.vec
	kept.list1=list1[which(list1 %in% celltype.order)]
	list2=celltypes.list$b.vec
	kept.list2=list2[which(list2 %in% celltype.order)]
	rowlabels=matrix(append(rep("dodgerblue4", length(kept.list1)), rep("deepskyblue", length(kept.list2))), nrow=1)
	
	med=(median(mat.to.plot, na.rm=TRUE))
	max=max(mat.to.plot, na.rm=TRUE)
	min=min(mat.to.plot, na.rm=TRUE)
	heatcols=gen_colors(min, med, max)

	if (length(ensembl.geneorder)==1)	{
		pdf(plotfile)
		x=t(mat.to.plot)
		print(x)
		par( mar = par( "mar" ) + c( 2, 4, 0, 0 ) )
		image(x, xaxt= "n", yaxt= "n", col=heatcols )
		celltype.annots=append(paste(kept.list1, "_", label1, sep=""), paste(kept.list2, "_", label2, sep=""))
		feat.order=paste(ensembl.geneorder, "_", mapgenes2chrstates, sep="")
		axis( 1, at=seq(0,1,length.out=nrow( x ) ), labels= feat.order, las= 2, cex.axis=.5)
		axis( 2, at=seq(0,1,length.out=ncol( x ) ), labels= celltype.annots, las= 2, cex.axis=.7 )
		dev.off()
		quit()	
	}

	## choose color scheme based on range of values
	allcolors=get_state_colors()
	
	## fill in colors corresponding to chromatin states
	default="mistyrose1"
	chromatin.labels=rep(default, ncol(mat.to.plot)) 
	for(currstate in unique(mapgenes2chrstates, na.rm=TRUE)) {
		if(!is.na(currstate)) {
			chromatin.labels[which(mapgenes2chrstates==currstate)]=allcolors[currstate]
		}
	}
	state.colors=FALSE
	if(state.colors) {
		columnlabels=matrix(chromatin.labels, ncol=1)
		colnames(columnlabels)=c("Chromatin_State")
	} else {
		columnlabels=matrix(nrow=0, ncol=0)
	}
	pdf(plotfile, width=14, height=14)
	columntext=rep("", ncol(mat.to.plot))
	if(useIds) {
		nspaces=ncol(mat.to.plot)/20
		tofill=seq(1,ncol(mat.to.plot), by=nspaces)
		columntext[tofill]=colnames(mat.to.plot)[tofill]
    }
	if(useNames) {
		currprefix=paste(metric, testtype, property, label1, label2, sep=".")
		common.genename.file=paste("all_pvals/", model, "/tables/", currprefix, "_matched.txt", sep="")
		nametable=read.table(common.genename.file)
		names(nametable)=c("GeneId", "State", "GeneName")
		columntext=as.vector(nametable$GeneName)
	}
	hm=heatmap.3(mat.to.plot, cexCol=0.5, RlabColor=get_ep_colors(rownames(mat.to.plot)), margins=c(9,12), keysize=.75, xlab=paste(ncol(mat.to.plot), "Genes"), ylab=paste(nrow(mat.to.plot), "Epigenomes"),  Colv=dend, Rowv=FALSE, dendrogram=dend.choice, scale="none", trace="none",RowSideColors=rowlabels, 
columnlabels, 
	col=heatcols, labCol=columntext, labRow=get_ep_names(rownames(mat.to.plot)))
	legend("topright", legend=c(label1, label2,"", "Multistate"), fill=c("dodgerblue4","deepskyblue", "white", default), border=FALSE, bty="n", y.intersp=0.7, cex=0.9)
	dev.off()
	savename=paste(plotdir, "exp.mat.to.plot", plottypesuffix, ".Rdata", sep="")
	print(paste("Saving mat.to.plot to ", savename, "...", sep=""))
	save(mat.to.plot, file=savename)
}
