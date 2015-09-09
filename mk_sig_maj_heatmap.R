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
source("funcs.R", chdir=T)
source("mkplots.R", chdir=T)
source("heatmap.3.R", chdir=T)

testing.bool=FALSE
if(testing.bool) {
	metric="perc"
	testtype="wilcox"
	correction="BY"
	property="type"
	group1="CellLine"
	group2="PrimaryCulture"
	label1="CellLine"
	label2="PrimaryCulture"
	model="core_gencode_v10_1kb"
	plottype="combinations"
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
matched=TRUE
useIds=FALSE
useNames=FALSE
suffix2=""
suffix=""
if(matched) {
	suffix="_matched"
}

set_variables(model, metadatafile,  genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label)
plottypesuffix=get_plottype_suffix(plottype)
list[reorder_by_assocstate, useDend, reorder_by_domstate, reorder_by_combinations] = get_plottype_bools(plottype)

dend=FALSE
dend.choice="none"

inputdir=paste("sig_majority/", model, sep="")
currprefix=paste(metric, testtype, correction, property, label1, label2, sep=".")
inputfile=paste(inputdir, "/", currprefix, "_matrix.txt", sep="")

plotdir=paste("plots/", model, "/", label1, ".", label2, "/", metric, "/", testtype, "/",correction, "/",  sep="")
dir.create(plotdir, recursive=TRUE, showWarnings=FALSE)
plotfile=paste(plotdir, "sig_maj_plot", suffix, suffix2, plottypesuffix, ".pdf", sep="")


#common.genename.file=paste("all_pvals/", model, "/gene_tables/", currprefix, "_matched.txt", sep="")
#print(paste("Checking for", common.genename.file, "..."))
#if(!file.exists(common.genename.file)) {
#	## make table file
#	dir.create(dirname(common.genename.file), recursive=TRUE, showWarnings=FALSE)
#	inputargs=paste(metric, testtype, correction, property, label1, label2, model)
#	setwd("all_pvals")
#	command=(paste("./mk_one_table.sh", inputargs))
#	print(command)
#	system(command)
#	setwd("..")
#}

mat=as.matrix(read.table(inputfile, header=TRUE, row.names=1))


sigpvals.order.file=paste(plotdir, "sig.pvals.rdata", sep="")


## for saving column order and associated chromatin states
geneorder.file=paste(plotdir, "sig_maj_geneorder", suffix, plottypesuffix,".txt", sep="")
chrstateorder.file=paste(plotdir, "sig_maj_chrstateorder", suffix, plottypesuffix,".txt", sep="")
celltypeorder.file=paste(plotdir, "sig_maj_celltypeorder", suffix, plottypesuffix,".txt", sep="")

library(gplots)
pair.label=get_pair_label(label1, label2)
celltypes.list=get_valid_celltypes(metric, property, pair.label)
kept.list1=celltypes.list$a.vec
kept.list2=celltypes.list$b.vec
celltype.order=c(kept.list1, kept.list2)

## get sig feats
sigpval.file=get_sigpval_file(model, metric, testtype, correction, property, label1, label2)
sigpvals=read.table(sigpval.file)
names(sigpvals)=c("feats", "pval")

## If no sig features, return. 
## If one sig features, special plotting necessary
if(nrow(sigpvals)==0) {
	print(paste("No significant features for ", featorder.file))
	return()
}

genes_in_order=get_genes_only(as.vector(sigpvals$feats))
names(genes_in_order)=as.vector(sigpvals$feats)

if(nrow(sigpvals)==1) {
	pdf(plotfile)
	x=t(mat[celltype.order, genes_in_order])
	## get coloring
	allcolors=get_state_colors()
	## choose color scheme based on range of values
	colors=allcolors[min(x, na.rm=TRUE):max(x, na.rm=TRUE)]

	if(min(x, na.rm=TRUE)==0) {
		## prepend white as many times as necessary to colors
		colors=append(rep("white", min(x[which(x!=0)])), colors)
	}
	par( mar = par( "mar" ) + c( 2, 4, 0, 0 ) )
	image(x, xaxt= "n", yaxt= "n", col=colors )
	celltype.annots=append(paste(kept.list1, "_", label1, sep=""), paste(kept.list2, "_", label2, sep=""))
	axis( 1, at=seq(0,1,length.out=nrow( x ) ), labels=sigpvals[1,1] , las= 2, cex.axis=.5)
	axis( 2, at=seq(0,1,length.out=ncol( x ) ), labels= celltype.annots, las= 2, cex.axis=.7 )
	dev.off()
	write.table(genes_in_order, file=geneorder.file, quote=FALSE, row.names=FALSE)
	chrstateorder=get_chrstates_only(sigpvals[1,1])
	write.table(chrstateorder, file=chrstateorder.file, quote=FALSE, row.names=FALSE)
	next()	
}

if(matched && useDend) {
	print("Ordering column ordering based on previous significant feat plots...")
	## load in all significant features
	featorder.file=paste(plotdir, "sigfeats_feats_featorder", suffix, ".txt", sep="")
	feat.order=as.vector(t(read.table(featorder.file, header=TRUE, colClasses="character")))
	## genes_in_order has multiple instances of some genes
	genes_in_order=get_genes_only(feat.order)
	col.names=unique(genes_in_order)
	#all_genes=unique(genes_in_order)
	#first_instances=match(all_genes, genes_in_order)
	#col.names=genes_in_order[first_instances]

	## use previous celltype ordering
	## load in celltype ordering
	celltypeorder.file=paste(plotdir, "sigfeats_feats_celltypeorder", suffix, ".txt", sep="")
	celltype.order=as.vector(t(read.table(celltypeorder.file, header=TRUE, colClasses="character")))
	row.names=celltype.order
} else {
	## reorder entire matrix based on dominant state
	currmat=mat[celltype.order,]
	
	## pick column ordering from all celltypes
	narm=TRUE
	print("Ordering columns from full heatmap...")
	## try to reorder columns
	col.order=get_order(currmat, 2, narm)	
	col.names=colnames(currmat)[col.order]
	
	## pick row ordering for each group separately
	print("Ordering rows...")
	mat1=mat[kept.list1, col.names]
	row1.ord=get_order(mat1, 1, narm)
	row1.names=rownames(mat1)[row1.ord]
		
	mat2=mat[kept.list2, col.names]
	row2.ord=get_order(mat2, 1, narm)
	row2.names=rownames(mat2)[row2.ord]
	
	row.names=append(row1.names, row2.names)
}

## combine into one matrix
ordered.mat=mat[row.names, col.names]

## map each gene to a chromatin state based on features
print("Mapping to colors...")
## get first index for each gene
first_inds=match(unique(genes_in_order), genes_in_order)
names(first_inds)=unique(genes_in_order)
sorted_inds=sort(first_inds, decreasing=FALSE)
## get features for these first occurrences of genes
feats=names(genes_in_order)[sorted_inds]
## get chrstates
mapgenes2chrstates=get_chrstates_only(feats)
names(mapgenes2chrstates)=get_genes_only(feats)
default="black"
fillvec=c("dodgerblue4","deepskyblue")
legendtext=c(label1, label2)

print("Reordering if necessary...")
if (reorder_by_assocstate) {
	new_order=order(as.numeric(mapgenes2chrstates))
} else if(reorder_by_domstate) {
	majdendfile=paste(plotdir, "sig_maj_plot_dend", plottypesuffix, ".rdata", sep="")
	dend=load_or_get_dend(majdendfile, ordered.mat, 2, TRUE,"dend")
	#save(list=c("ordered.mat"), file="tmp.mat.1.rdata")
	dend.choice="column"
	new_order=1:ncol(ordered.mat)
} else if(reorder_by_combinations) {
	load(sigpvals.order.file)
	new_order=sigpvals.order
} else {
	new_order=1:ncol(ordered.mat)
}
ordered.mat=ordered.mat[,new_order]
mapgenes2chrstates=mapgenes2chrstates[new_order]

print("Generating final heatmap...")
## make row colors that show relevant grouping
rowlabels=matrix(append(rep("dodgerblue4", length(kept.list1)), rep("deepskyblue", length(kept.list2))), nrow=1)
rownames(rowlabels)=c("Group")

colors=get_chrstate_colors(min(ordered.mat, na.rm=TRUE):max(ordered.mat, na.rm=TRUE))

state.colors=FALSE
if(state.colors) {
	chromatin.labels=get_chrstate_colors(mapgenes2chrstates)
	clab=matrix(chromatin.labels, ncol=1)
	colnames(clab)=c("Chromatin_state")
} else {
	clab=matrix(nrow=0, ncol=0)
}

columntext=rep("", ncol(ordered.mat))
if(useIds){
	nspaces=ncol(ordered.mat)/40
	tofill=seq(1,ncol(ordered.mat), by=nspaces)
	columntext[tofill]=colnames(ordered.mat)[tofill]
}
if(useNames) {
	columntext=as.vector(get_gene_symbols(colnames(ordered.mat)))
}
print(paste("Plotting to", plotfile, "..."))
pdf(plotfile, width=14, height=14)
hm=heatmap.3(ordered.mat, cexCol=.5, RlabColor=get_ep_colors(rownames(ordered.mat)), margins=c(9, 12), keysize=.75, ylab=paste(nrow(ordered.mat), "Epigenomes"), xlab=paste(ncol(ordered.mat), "genes"), 
	Colv=dend, Rowv=FALSE, dendrogram=dend.choice, scale="none", trace="none",RowSideColors=rowlabels, 
	ColSideColors=clab, 
	col=colors, labRow=get_ep_names(rownames(ordered.mat)), labCol=columntext)
legend("topright", legend=legendtext, fill=fillvec, border=FALSE, bty="n", y.intersp=0.7, cex=0.9)
dev.off()

print(paste("Writing to", geneorder.file))
write.table(rownames(ordered.mat), file=celltypeorder.file, quote=FALSE, row.names=FALSE)
write.table(colnames(ordered.mat), file=geneorder.file, quote=FALSE, row.names=FALSE)
write.table(mapgenes2chrstates, file=chrstateorder.file, quote=FALSE, row.names=FALSE)

save(ordered.mat, file=paste(plotdir, "sig_maj_ordered.mat", plottypesuffix, ".rdata", sep=""))

