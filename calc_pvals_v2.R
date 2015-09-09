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

source("funcs.R", chdir=T)
source("setvars.R", chdir=T)
source("clust_funcs.R", chdir=T)
library("ggplot2")
## calculates p-vals for expression differences overall and in specific sublabels.str
## Based on Mann-Whitney p-val
testone=FALSE
firstrun=FALSE
if(testone) {
    property="specialgi"
    a.option="ESC"
    b.option="GI"
    test="wilcox"
    correction="fdr"
	generegions_label="gencode_v10_1kb"
    model=paste0("core_", generegions_label)
	firstrun=FALSE
	metadatafile="data/final_celltype_metadata.txt"
	genefile="data/gencode_genes_1kb_up.txt"
	covariate_mat_file="data/cov.mat.txt"
	map_covariates_file="data/map_vars_covariates.txt"
	expfile="data/57epigenomes.RPKM.pc"
	state_annotations_file="data/core_annotation.txt"
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
## vars needed to be set
set_variables(model, metadatafile,  genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label)
a.label=a.option
b.label=b.option
metric="perc"
pair=get_pair_label(a.label, b.label)
plottype="domstate"

suffix="_matched"
##genelabels.str should be read in from annot_genegroups/${model}/ dir

resultsdir="exp.pvals/"
dir.create(resultsdir, showWarnings=FALSE, recursive=TRUE)

resid.result.table=c()
orig.result.table=c()

plottypesuffix=get_plottype_suffix(plottype)
outfile=paste(resultsdir, model, "_", plottype, "_percent_all.txt", sep="")
resid.tablefile=paste(resultsdir, model, "_", plottype, "_percent_all_resid_table.txt", sep="")
orig.tablefile=paste(resultsdir, model, "_", plottype, "_percent_all_orig_table.txt", sep="")
clust.outfile=paste(resultsdir, model, "_", plottype, "_clusters_percent_all.txt", sep="")

if(!file.exists(resid.tablefile) && !file.exists(orig.tablefile)) {
	firstrun=TRUE
}
if(firstrun) {mode="w"} else {mode="a"}
file.conn=file(outfile, mode)
clust.file.conn=file(clust.outfile, mode)
resid.table.file.conn=file(resid.tablefile, mode)
orig.table.file.conn=file(orig.tablefile, mode)

get_expression_qvals=function(expmat, dim, celltypes1, celltypes2, method="BH") {
	results=apply(expmat, dim, function(expvals) {
			if(all(is.na(expvals))) { return(NA)} else {wilcox.test(rm.nas(as.vector(expvals[celltypes1])), rm.nas(as.vector(expvals[celltypes2])), alternative="two.sided")$p.val }
			})
	q.vals=p.adjust(rm.nas(results), method="BH")
	return(q.vals)
}

plotdir=get_full_plotdir(model, a.label, b.label, test, metric, correction)


## load genes in order
geneorderfile=get_geneorder_file(plotdir, suffix, plottypesuffix)
if(!file.exists(geneorderfile)) {
    stop(paste0(geneorderfile, " does not exist. Please run earlier parts of ChromDiff analysis before running calc_pvals.R"))
}
genes=as.vector(t(read.table(geneorderfile, skip=1)))
## drop end of gene id
ensemblids=get_ensembl.ids(genes)
residsdir=paste0(resultsdir, model, "/resids/")
dir.create(residsdir, showWarnings=FALSE, recursive=TRUE)
expdatafile=get_expdata_file(plotdir, plottypesuffix)
allpvalfile=get_pval_file(model, metric, test, correction, property, a.label, b.label)
## CALCULATE EXPRESSION DIFFERENCES ON ALL SIG FEATS
if(!file.exists(expdatafile) || !file.exists(allpvalfile)) {
	stop(paste0("Either ", expdatafile, " or ", allpvalfile, " does not exist."))
}
resid.currcol=c()
orig.currcol=c()

writeLines(paste("Results for ", plotdir, "...", sep=""), file.conn)

pval.table=read.table(allpvalfile)
names(pval.table)=c("featnames", "corr.pvals")

## figure out which genes were associated with sig features
keep.inds=which(pval.table[,"corr.pvals"]<.05)
keep.feats=pval.table[keep.inds,"featnames"]
gencode.genes=unique(get_genes_only(as.vector(keep.feats)))
ensembl.genes=unique(get_ensembl.ids(gencode.genes))

## get groupings of epigenomes
celltypes.list=get_valid_celltypes(metric, property, pair)
kept.list1=celltypes.list$a.vec
kept.list2=celltypes.list$b.vec

## resids.mat and logexp.mat should have same number of rows (all eps with exp data)
fullresids.mat=get_exp_resids(property, residsdir)
## fill in mats for needed eps
resids.mat=fill_nas(fullresids.mat, c(kept.list1, kept.list2), ensembl.genes)
logexp.mat=get_logexp(c(kept.list1, kept.list2), ensembl.genes)
nonsig.genes=get_other_genes(ensembl.genes)
nonsig.mat=get_logexp(c(kept.list1, kept.list2), nonsig.genes)
nonsig.resids=fill_nas(fullresids.mat, c(kept.list1, kept.list2), nonsig.genes)

if(is.null(logexp.mat)) {
	stop("No expression data available for specified epigenomes.") 
}
## test resids for expression differences among sig genes and nonsig genes. print results.
q.vals=get_expression_qvals(resids.mat, 2, kept.list1, kept.list2)
writeLines(paste(length(which(q.vals<.05)), "out of", length(q.vals), "sig genes (", length(which(q.vals<.05))*100/length(q.vals), "%) show significant expression based on covariate corrected values"), file.conn)
resid.currcol=append(resid.currcol, c(length(which(q.vals<.05)), length(q.vals)))
q.vals=get_expression_qvals(nonsig.resids, 2, kept.list1, kept.list2)
writeLines(paste(length(which(q.vals<.05)), "out of", length(q.vals), "non-sig genes (", length(which(q.vals<.05))*100/length(q.vals), "%) show significant expression based on original values"), file.conn)
resid.currcol=append(resid.currcol, c(length(which(q.vals<.05)), length(q.vals)))

## process results into table 
resid.m=matrix(resid.currcol, nrow=1)
colnames(resid.m)=c("sigexp.siggenes", "numsiggenes", "sigexp.nonsiggenes", "numnonsiggenes")
rownames(resid.m)=paste0(pair, ".", correction)
resid.result.table=rbind(resid.result.table, resid.m)

## test original uncorrected expression values for differences in sig genes and nonsig genes. print results
q.vals=get_expression_qvals(logexp.mat, 2, kept.list1, kept.list2)
writeLines(paste(length(which(q.vals<.05)), "out of", length(q.vals), "sig genes (", length(which(q.vals<.05))*100/length(q.vals), "%) show significant expression based on original values"), file.conn)
orig.currcol=append(orig.currcol, c(length(which(q.vals<.05)), length(q.vals)))
q.vals=get_expression_qvals(nonsig.mat, 2, kept.list1, kept.list2)
writeLines(paste(length(which(q.vals<.05)), "out of", length(q.vals), "non-sig genes (", length(which(q.vals<.05))*100/length(q.vals), "%) show significant expression based on original values"), file.conn)
orig.currcol=append(orig.currcol, c(length(which(q.vals<.05)), length(q.vals)))	

## process results into table
orig.m=matrix(orig.currcol, nrow=1)
colnames(orig.m)=c("sigexp.siggenes", "numsiggenes", "sigexp.nonsiggenes", "numnonsiggenes")
rownames(orig.m)=paste0(pair, ".", correction)
orig.result.table=rbind(orig.result.table, orig.m)

clusters=TRUE
if(clusters) {
	##repeat with specific gene labels.str
	vert.ints.file=get_vert.ints.file(plotdir, plottypesuffix)
    if(!file.exists(vert.ints.file)) {
    	stop(paste0(vert.ints.file, " does not exist."))
    }
	# load clusterings
    load(vert.ints.file)
    if(!is.list(vert.ints.noreps)) {
        stop("vert.ints.noreps is not a list.")
    }
    vert.ints.noreps=order.list.by.first(vert.ints.noreps)

	## split exp vals into groupings
    clustcount=1
	writeLines(paste("Results for ", plotdir, "...", sep=""), clust.file.conn)
    for(clustind in 1:length(vert.ints.noreps)) {
        currclust.inds=vert.ints.noreps[[clustind]]
        if((currclust.inds[2]-currclust.inds[1])<=minclustsize) {
            next
        }
        print(paste("Calculating clust", clustcount, "..."))
		## check what percent of this group's genes have sig differences
        genegroup=unique(ensemblids[currclust.inds[1]:currclust.inds[2]])
		genegroup=rm.nas(as.vector(genegroup))
		currclust.logexp.mat=logexp.mat[, genegroup]
		q.vals=get_expression_qvals(currclust.logexp.mat, 2, kept.list1, kept.list2)
		writeLines(paste(length(which(q.vals<.05)), "out of", length(q.vals), "in gene subgroup", clustcount,  "(", length(which(q.vals<.05))*100/length(q.vals), "%) show significant expression based on original values"), clust.file.conn)

		## calculate overal difference from original expression
		group1.vec=rm.nas(as.vector(logexp.mat[kept.list1, genegroup]))
		group2.vec=rm.nas(as.vector(logexp.mat[kept.list2, genegroup]))
		if(length(group1.vec)==0 || length(group2.vec)==0) {
			writeLines(paste("Skipping gene group", clustcount, "because no gene expression data..."), clust.file.conn)
		} else {
			orig.pval= wilcox.test(group1.vec, group2.vec, alternative="two.sided")$p.val
			writeLines(paste("P-value of uncorrected expression for gene group", clustcount, "is", orig.pval), clust.file.conn)
		}

		## calculate overal difference from residual expression
		resids.group1.vec=rm.nas(as.vector(resids.mat[kept.list1, genegroup]))
		resids.group2.vec=rm.nas(as.vector(resids.mat[kept.list2, genegroup]))
		if(length(resids.group1.vec)==0 || length(resids.group2.vec)==0) {
			writeLines(paste("Skipping gene group", clustcount, "because no gene expression data..."), clust.file.conn)
		} else {
			## write pval info
			resids.pval=wilcox.test(resids.group1.vec, resids.group2.vec, alternative="two.sided")$p.val
			writeLines(paste("P-value of covariate corrected expression for gene group", clustcount, "is", resids.pval), clust.file.conn)
		}
		clustcount=clustcount+1
	}
} 
writeLines("", clust.file.conn)
writeLines("", file.conn)
write.table(resid.result.table, resid.table.file.conn, col.names=firstrun, quote=FALSE, sep="\t")
write.table(orig.result.table, orig.table.file.conn, col.names=firstrun, quote=FALSE, sep="\t")
close(orig.table.file.conn)
close(resid.table.file.conn)
close(clust.file.conn)
close(file.conn)

