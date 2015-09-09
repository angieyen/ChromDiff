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
source("clust_funcs.R", chdir=T)
test=FALSE
sigfeat_plot=FALSE
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
	heightcutoff=80
	plottype.str="domstate"
    metadatafile="data/final_celltype_metadata.txt"
    genefile="data/gencode_genes_1kb_up.txt"
    covariate_mat_file="data/cov.mat.txt"
    map_covariates_file="data/map_vars_covariates.txt"
    expfile="data/57epigenomes.RPKM.pc"
    state_annotations_file="data/core_annotation.txt"
	generegions_label="gencode_v10"
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
	heightcutoff=as.numeric(args[10])
	## options are assocstate usedend and domstate
	plottype.str=args[11]
    metadatafile=args[12]
    genefile=args[13]
    covariate_mat_file=args[14]
    map_covariates_file=args[15]
    expfile=args[16]
    state_annotations_file=args[17]
	generegions_label=args[18]
}
set_variables(model, metadatafile,  genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label)


### get parameters used by all/multiple plots
plotdir=get_full_plotdir(model, label1, label2, testtype, metric, correction)

suffix="_matched"

plottypesuffix=get_plottype_suffix(plottype.str)
list[reorder_by_assocstate, reorder_by_dend, reorder_by_domstate, reorder_by_combinations] = get_plottype_bools(plottype.str)

majdendfile=paste(plotdir, "sig_maj_plot_dend", plottypesuffix, ".rdata", sep="")
chrstateorder.file=paste(plotdir, "sig_maj_chrstateorder", suffix, plottypesuffix,".txt", sep="")

##output clustered plots at the following files
plot.file=paste(plotdir,"sigfeats_feats_clust", suffix, plottypesuffix, ".pdf", sep="")
exp.plotfile=paste(plotdir, "sig_exp_plot_clust", suffix, plottypesuffix, ".pdf", sep="")
majstate.plotfile=paste(plotdir, "sig_maj_plot_clust", suffix, plottypesuffix, ".pdf", sep="")

linewidth=5

## load in celltype ordering
pair.label=get_pair_label(label1, label2)
celltypes.list=get_valid_celltypes(metric, property, pair.label)
list1=celltypes.list$a.vec
list2=celltypes.list$b.vec

## make row colors that show relevant grouping
rowlabels=matrix(append(rep("dodgerblue4", length(list1)), rep("deepskyblue", length(list2))), nrow=1)
row.side.labels=c(paste0(label1, " (", length(list1), " epigenomes)"), paste0(label2, " (", length(list2), " epigenomes)"))

## write text for group 1 and group 2
## calc positions of group labels
total=length(list1)+length(list2)
pos1=total-length(list1)/2
pos2=length(list2)/2
row.side.fracs=c(pos1/total, pos2/total)

cexRowSideLabels=c(2.5,2.5)
cexColSideLabels=c(2)
cexlab=3
cex.row=.8
rightmargin=3
margins=c(4, rightmargin)

## set position parameters so key is on top right and margins/size divisions make sense
lmat=rbind(c(0, 3, 4), c(2, 1, 0))
lwid=c(.2, 6, 1.5)
lhei=c(0.75,7)
horiz.ints=c(length(list2)+.5)
plotwidth=25
plotheight=13
	

################# MAJORITY STATE PLOT ################
maj.rdata=paste(plotdir, "sig_maj_ordered.mat", plottypesuffix, ".rdata", sep="")
load(maj.rdata)
dend=FALSE
dend.choice="none"

roworder=rownames(ordered.mat)

if(reorder_by_dend || reorder_by_domstate || reorder_by_combinations) {
	minclustfrac=.05
	minclustsize=minclustfrac*ncol(ordered.mat)
} else if (reorder_by_assocstate) {
	minclustsize=1
}

ls.results=get_vert.ints(plottype.str, plotdir, ordered.mat)
vert.ints=ls.results$vert.ints
ordered.mat=ls.results$mat.to.plot

## load in chromatin state annotations 
chrstateorder.file=paste(plotdir, "sig_maj_chrstateorder", suffix, plottypesuffix, ".txt", sep="")
chrstates_norep_order=as.vector(t(read.table(chrstateorder.file, header=TRUE)))
## to order by associated chr state
if(reorder_by_assocstate){
	## find the first instance of each chromatin state
	first.inst=sapply(min(chrstates_norep_order):max(chrstates_norep_order), function(x) {which(chrstates_norep_order==x)[1]})
	first.inst=first.inst[which(!is.na(first.inst))]
	vert.ints.noreps=list()
	for(ind in 1:length(first.inst)) {
		if(ind<length(first.inst)) {
			vert.ints.noreps[[ind]]=c(first.inst[ind], first.inst[ind+1]-1)
		} else {
			vert.ints.noreps[[ind]]=c(first.inst[ind], length(chrstates_norep_order))
		}
	}
} else if(reorder_by_domstate) {
	## order by majority state (pick one feature for each gene)
	## feat plot already didn't have features with repeating genes, so the intercepts stay the same
	vert.ints.noreps=vert.ints
	dend=load_or_get_dend(majdendfile, ordered.mat, 2, TRUE, "dend")
	hclust.file=paste(plotdir, "hclust", plottypesuffix, ".pdf", sep="")
	dend.choice="column"
	pdf(hclust.file)
	plot(dend)
	dev.off()
} else if(reorder_by_combinations) {
	vert.ints.noreps=vert.ints
}

## choose color scheme based on range of values
allcolors=get_state_colors()
colors=allcolors[min(ordered.mat, na.rm=TRUE):max(ordered.mat, na.rm=TRUE)]
pdf(majstate.plotfile, width=plotwidth, height=plotheight)
## write label for Chromatin State Coloring
plot.title="Dominant chromatin state"
print(paste("Making", plot.title, "plot..."))
par(cex.main=3, cex.lab=cexlab)
columntext=rep("", ncol(ordered.mat))
## set x label axis for gene plots (expression and majority state)
xlab=paste(ncol(ordered.mat), "significant genes")
#save(list=ls(all=TRUE), file="sex.rdata")
hm=heatmap.3(ordered.mat, main=plot.title,  margins=margins, key=FALSE,  cexRow=cex.row, col=colors, xlab=xlab, 
	#ylab=ylab, 
	Colv=dend, Rowv=FALSE, dendrogram=dend.choice, scale="none", trace="none",
	row.side.height.fraction=0.1,
	col.side.height.fraction=0.1,
	labRow=get_ep_names(rownames(ordered.mat)), RlabColor=get_ep_colors(rownames(ordered.mat)), RowSideColors=rowlabels, RowSideLabels=row.side.labels, RowSideFracs=row.side.fracs, cexRowSideLabels=cexRowSideLabels,
	labCol=columntext, 
	ColSideFracs=col.side.fracs, cexColSideLabels=cexColSideLabels,
	lmat=lmat, lwid=lwid, lhei=lhei, 
	add.expr=final.drawing(horiz.ints, vert.ints.noreps, linewidth, minclustsize, nrow(ordered.mat)))
a=dev.off()

vert.ints.file=get_vert.ints.file(plotdir, plottypesuffix)
print(paste("Writing to", vert.ints.file))
save(list=c("vert.ints.noreps", "minclustsize"), file=vert.ints.file)

################## GENE EXPRESSION PLOT ###################################
plot.title="Gene expression differences"

dend=FALSE
dend.choice="none"

## now load in expression data to plot those with boxes/lines
exp.rdata=paste(plotdir, "exp.mat.to.plot", plottypesuffix, ".Rdata", sep="")
if(file.exists(exp.rdata)) {
	print(paste("Making", plot.title, "plot..."))
	load(exp.rdata)
	
	exp.mat.to.plot=mat.to.plot[roworder,]
	
	## find relevant cluster vertical intercepts for gene expression
	if(reorder_by_dend) {
		ls=get_vert.ints_dend(plotdir, exp.mat.to.plot, noreps=TRUE)
		vert.ints.noreps=ls$vert.ints.noreps
	} else if(reorder_by_assocstate){
		first.inst=sapply(min(chrstates_norep_order):max(chrstates_norep_order), function(x) {which(chrstates_norep_order==x)[1]})
		first.inst=first.inst[which(!is.na(first.inst))]
		vert.ints.noreps=list()
		for(ind in 1:length(first.inst)) {
			if(ind<length(first.inst)) {
				vert.ints.noreps[[ind]]=c(first.inst[ind], first.inst[ind+1]-1)
			} else {
				vert.ints.noreps[[ind]]=c(first.inst[ind], length(chrstates_norep_order))
			}
		}
	} else if(reorder_by_domstate || reorder_by_combinations) {
		## feat plot already didn't have features with repeating genes, so the intercepts stay the same
		vert.ints.noreps=vert.ints
	}

	if(reorder_by_domstate) {
		dend=load_or_get_dend(majdendfile, exp.mat.to.plot, 2, TRUE, "dend")
		dend.choice="column"
	}
	#print(vert.ints.noreps)
	

	## set chromatin state coloring vector
	col.side.fracs=c(NA)

	## make colorscale for exp mat
	med=(median(exp.mat.to.plot, na.rm=TRUE))
	max=max(exp.mat.to.plot, na.rm=TRUE)
	min=min(exp.mat.to.plot, na.rm=TRUE)
	heatcols=gen_colors(min, med, max)

	pdf(exp.plotfile, width=plotwidth, height=plotheight)
	## make expression plot
	par(cex.main=3, cex.lab=cexlab)
	rightmargin=100/nrow(exp.mat.to.plot)
	if(reorder_by_dend || reorder_by_domstate || reorder_by_combinations) {
		minclustfrac=.05
		minclustsize=minclustfrac*ncol(exp.mat.to.plot)
	} else if (reorder_by_assocstate) {
		minclustsize=1
	}
	KeyValueName="Gene expression"
	hm=heatmap.3(exp.mat.to.plot, main=plot.title, cexRow=cex.row,  margins=margins,  xlab=xlab, 
		#ylab=ylab, 
		col=heatcols, keysize=0.5,  
		Colv=dend, Rowv=FALSE, dendrogram=dend.choice, scale="none", trace="none",
		row.side.height.fraction=0.1,
		col.side.height.fraction=0.1,
		RowSideColors=rowlabels,  RlabColor=get_ep_colors(rownames(exp.mat.to.plot)), RowSideLabels=row.side.labels, RowSideFracs=row.side.fracs, labRow=get_ep_names(rownames(exp.mat.to.plot)), cexRowSideLabels=cexRowSideLabels,
		labCol=columntext,ColSideFracs=col.side.fracs,cexColSideLabels=cexColSideLabels,
		KeyValueName=KeyValueName,
		lmat=lmat, lwid=lwid, lhei=lhei,
		add.expr=final.drawing(horiz.ints, vert.ints.noreps, linewidth, minclustsize, nrow(exp.mat.to.plot)))
	a=dev.off()
} else {
	print(paste("Expression data not available for", exp.plotfile))
}


##### Significant features plot #####
if(sigfeat_plot) {
	mat.file=paste(plotdir,"sigfeats_feats_mat.to.plot_matched.rdata", sep="")

	## matrix saved as mat.to.plot
	load(mat.file)
	feat.mat.to.plot=mat.to.plot[roworder,]
	total_cols=ncol(feat.mat.to.plot)
	total_rows=nrow(feat.mat.to.plot)
	if(total_cols==1) {
		stop("No clustered plots due to only one sig feature");
	}
	if (reorder_by_domstate || reorder_by_combinations) {
		xlab=paste(ncol(feat.mat.to.plot), "significant genes")
	} else {
		xlab=paste(ncol(feat.mat.to.plot), "significant features")
	}
	print("Loading in celltype, gene, and chromatin state ordering...")


	## fill in colors corresponding to chromatin states
	colstates = as.character(get_chrstates_only(colnames(feat.mat.to.plot)))
	chromatin.labels=get_chrstate_colors(colstates)

	print("Making gene group labels...")
	genegroup.legends=c("Group 1", "Group 1 matched", "Unique")
	genegroup.colors=c("mistyrose1", "mistyrose3", "mistyrose4")
	if(reorder_by_assocstate) {
		clab=matrix(chromatin.labels, ncol=1)
		colnames(clab)=c("Chromatin state")
	} else {
		clab=matrix(nrow=0, ncol=0)
	}
	all.col.legends=c(label1, label2, "",genegroup.legends)
	group.cols=c("dodgerblue4", "deepskyblue")
	all.col.fill=c(group.cols, "white", genegroup.colors)


	## fill in column labels appropriately
	useIds=FALSE
	useNames=FALSE
	if(useIds) {
		nspaces=ncol(feat.mat.to.plot)/20
		tofill=seq(1,ncol(feat.mat.to.plot), by=nspaces)
		columntext[tofill]=colnames(feat.mat.to.plot)[tofill]
	}
	if(useNames) {
		currprefix=paste(metric, testtype, property, label1, label2, sep=".")
		common.genename.file=paste("all_pvals/", model, "/tables/", currprefix, "_matched.txt", sep="")
		nametable=read.table(common.genename.file)
		names(nametable)=c("GeneId", "State", "GeneName")
		columntext=as.vector(nametable$GeneName)
	}

	heatcols=colorRampPalette(c("blue", "white", "red"))(100)
	# make x label axis as # features for sig feats plots
	# make y label axis for number of epigenomes
	#ylab=paste(nrow(feat.mat.to.plot), "epigenomes")
	rlabcolors=get_ep_colors(rownames(feat.mat.to.plot))


	pdf(plot.file, width=plotwidth, height=plotheight)
	plot.title="Relative enrichment for significant features"
	print(paste("Making", plot.title, "plot..."))
	par(cex.main=3, cex.lab=cexlab)
	if(reorder_by_dend || reorder_by_domstate || reorder_by_combinations) {
		minclustfrac=.05
		minclustsize=minclustfrac*ncol(feat.mat.to.plot)
	} else if (reorder_by_assocstate) {
		minclustsize=1
	}

	hm=heatmap.3(feat.mat.to.plot, main=plot.title, cexRow=cex.row, margins=margins,  xlab=xlab, 
		#ylab=ylab,
		keysize=0.5,
		Colv=FALSE, Rowv=FALSE, dendrogram="none", scale="none", trace="none",
		col.side.height.fraction=0.1,
		row.side.height.fraction=0.1,
		labRow=get_ep_names(rownames(feat.mat.to.plot)),
		RlabColor=rlabcolors, RowSideColors=rowlabels, RowSideLabels=row.side.labels, RowSideFracs=row.side.fracs, labRow=epnames, cexRowSideLabels=cexRowSideLabels,    
		ColSideColors=clab, 
		col=heatcols, labCol=columntext, cexColSideLabels=cexColSideLabels,
		lmat=lmat, lwid=lwid, lhei=lhei,
		add.expr=final.drawing(horiz.ints, vert.ints, linewidth, minclustsize, total_rows))
	a=dev.off()
}
