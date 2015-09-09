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
source("heatmap.3.R", chdir=T)
suppressMessages(library(gplots, quietly=TRUE))


make.plots= function(property, group.a, group.b, a.label, b.label,  metric, test, test.correction, model, no_covariate_correction, sigfeat.plots=FALSE, figure1=FALSE, pval.plots=TRUE, random.plots=FALSE) {
   	print(test.correction) 
	save(list=ls(all=TRUE), file="tmp.rdata")
	sigpval.file=get_sigpval_file(model, metric, test, test.correction, property, a.label, b.label)
	if(!file.exists(sigpval.file)) {
		return()
	}

	## if we have significant features, make plots for them
	sig.feats=as.vector(rownames(read.table(sigpval.file, row.names=1)))
	allpval.file=get_pval_file(model, metric, test, test.correction, property, a.label, b.label)
	corrpval.mat=read.table(allpval.file, row.names=1)
	corrected.pvals=as.vector(corrpval.mat[,1])
	names(corrected.pvals)=rownames(corrpval.mat)
    
    resid.mat=get_residuals(model, metric, property, no_covariate_correction)
    valid_celltypes_a=group.a[which(group.a %in% colnames(resid.mat))]
    valid_celltypes_b=group.b[which(group.b %in% colnames(resid.mat))]
    matrix.a=resid.mat[,valid_celltypes_a]
    matrix.b=resid.mat[,valid_celltypes_b]
	
	plotdir=get_plotdir(model, a.label, b.label)
	subdir=get_plot_subdir(test, metric, test.correction)
	statevec=get_all_chrstates()

	print("Starting to make plots...")

	if(figure1) {
    	make.subsetfeat.heatmaps(matrix.a, matrix.b, a.label, b.label, sig.feats,  corrected.pvals, plotdir, subdir, statevec, FALSE, model)
	}

	if(random.plots){
    	limit=1000
    	top.feats=names(corrected.pvals)[1:limit]
    	make.sigfeat.heatmaps(matrix.a, matrix.b, a.label, b.label, top.feats, plotdir, subdir, TRUE, model, filebase="top1000feats")               
	}

	if(pval.plots) {
    	mk.all.pval.plots(corrected.pvals[sig.feats], plotdir, subdir, model)
	}

	## usually only plot matched heatmaps
	if(sigfeat.plots) {
    	matched=TRUE
    	make.sigfeat.heatmaps(matrix.a, matrix.b, a.label, b.label, sig.feats, plotdir, subdir, matched, model, filebase="sigfeats")
	}
	genesize.plots=TRUE
	if(genesize.plots) {
		make.genesize.plots(sig.feats, plotdir, subdir)
	}


}
make.genesize.plots <- function(sig.feats, plotdir, subdir) {
	# make density plots for genesizes of sig, all, and nonsig genes	
	plot.subdir=paste(plotdir, subdir, sep="")
	dir.create(plot.subdir, recursive=TRUE, showWarnings=FALSE)	
	
	sigplotfile=paste0(plot.subdir, "genesize_sig_dist.pdf")
	allplotfile=paste0(plot.subdir, "genesize_all_dist.pdf")
	nonsigplotfile=paste0(plot.subdir, "genesize_nonsig_dist.pdf")
	
	sig_geneids=unique(get_genes_only(sig.feats))
	sig_genesizes=get_genesizes(sig_geneids)
	#plot.genesizes(sig_genesizes, sigplotfile)

	geneids=get_all_genes()
	genesizes=get_genesizes(geneids)
	#plot.genesizes(genesizes, allplotfile)
	
	nonsig_geneids=geneids[which(!(geneids %in% sig_geneids))]
	nonsig_genesizes=get_genesizes(nonsig_geneids)
	#plot.genesizes(nonsig_genesizes, nonsigplotfile)

	result=wilcox.test(sig_genesizes, nonsig_genesizes, paired=FALSE)
	plot.comp.genesizes(sig_genesizes, nonsig_genesizes, plot.subdir, result$p.value)

}

plot.genesizes=function(genesizes, file) {
	library("ggplot2")
	dat=data.frame(size=genesizes)
	ggplot(dat, aes(x=size)) + geom_density()
	ggsave(file)
}

plot.comp.genesizes=function(sig_genesizes, nonsig_genesizes, plot.subdir, pval) {
	library("ggplot2")
	allsizes=c(sig_genesizes, nonsig_genesizes)
	labels=c(rep("sig", times=length(sig_genesizes)), rep("nonsig", times=length(nonsig_genesizes)))
	dat=data.frame(size=allsizes, label=labels)
	#ggplot(dat, aes(x=size, colour=label))+geom_density()
	p0=ggplot(dat, aes(y=size, x=label))+geom_boxplot()+guides(fill=FALSE) + ggtitle(paste("P-value is", pval))
	file=paste0(plot.subdir, "genesize_comp.pdf")
	ggsave(file)
	ylim1=boxplot.stats(dat$size)$stats[c(1,5)]
	p1=p0+coord_cartesian(ylim=ylim1*1.25)
	file=paste0(plot.subdir, "genesize_comp_zoom.pdf")
	ggsave(file)
}
plot.top.n <- function(orig.mat, ranks, limit, reorderRows, reorderCols, plotdir, filebase, heatcols, model, colorStates) {
	save(list=ls(all=TRUE), file="tmp.brain.gi.rdata")
	print(paste("Plotting top", limit, "features..."))
	## for any values that are not in top n, set value to NA
	default=NA
	comp_func <- is.na
	indices.to.keep=which(ranks<limit)
	top.vals=rep(default, length(ranks))
	top.vals[indices.to.keep]=orig.mat[indices.to.keep]
	top.mat=matrix(top.vals, nrow=dim(orig.mat)[1], dimnames=list(rownames(orig.mat), colnames(orig.mat)))
	
	too.little.data.msg=paste("Not enough data in the matrix to plot", filebase)
	
	## identify cols and rows you are throwing out because there are no values in the top n
	if(ncol(top.mat)>1) {
		keepcols = which(apply(top.mat, 2, function(col) { return(length(which(comp_func(col)))!=length(col))}))
		if(length(keepcols)==0) { print(too.little.data.msg); return()}
		filtered.mat=top.mat[, keepcols]
		if (nrow(top.mat)==1 || length(keepcols)==1) {
			filtered.mat=matrix(filtered.mat, nrow=nrow(top.mat), dimnames=list(rownames(top.mat), colnames(top.mat)[keepcols]))
		}
	} else {
		filtered.mat=top.mat
	}
	if(nrow(filtered.mat)>1) {
		keeprows= which(apply(filtered.mat, 1, function(row) { return(length(which(comp_func(row)))!=length(row))}))
		if(length(keeprows)==0) { print(too.little.data.msg); return()}
		final.mat=filtered.mat[keeprows,]
		if (ncol(filtered.mat)==1 || length(keeprows)==1) {
			final.mat=matrix(final.mat, ncol=ncol(filtered.mat),  dimnames=list(rownames(filtered.mat)[keeprows], colnames(filtered.mat)))
		}
	} else {
		final.mat=filtered.mat
	}
	
	outfile=paste(plotdir, filebase, ".pdf", sep="")
	## put back in original values
	final.filled.mat=orig.mat[rownames(final.mat), colnames(final.mat)]
	## put in "p-values" of 1 if NA
	final.filled.mat[which(is.na(final.filled.mat))]=1
	if(!is.matrix(final.filled.mat)) { final.filled.mat=matrix(final.filled.mat, nrow=nrow(final.mat), dimnames=list(rownames(final.mat), colnames(final.mat))) }
	rownames(final.filled.mat) = mapStateNumToName(rownames(final.filled.mat))
	sigpvals.order=c()
	if(nrow(final.filled.mat)==1 || ncol(final.filled.mat)==1) { 
		plot_vector_as_hm(outfile, final.filled.mat,rownames(final.filled.mat), colnames(final.filled.mat))
		sigpvals.order=colnames(final.filled.mat)
	} else if(is.matrix(final.filled.mat)) {
		if(reorderCols && reorderRows) {
			dend="both"
		} else if (reorderCols) {
			dend="column"
		} else if (reorderRows) {
			dend="row"
		} else {
			dend="none"
		}

		pdf(outfile)
		hmlist = tryCatch({
			hm=heatmap.2(final.filled.mat, density.info="none", na.rm=TRUE, na.color="white", trace="none", cexCol=0.1, dendrogram=dend, col=heatcols, Rowv=reorderRows, Colv=reorderCols, scale="none")

			if(colorStates) {
				dev.off()
				reord.mat=final.filled.mat[hm$rowInd, hm$colInd]

				reord.mat[which(reord.mat==1)]=NA
				stateNames=rownames(reord.mat)
				for(stateName in stateNames) {
					colinds=which(!(is.na(reord.mat[stateName,])))
					reord.mat[stateName, colinds] = as.numeric(mapStateNameToNum(stateName))
				}
				minstatenum=min(as.numeric(mapStateNameToNum(stateNames)))
				maxstatenum=max(as.numeric(mapStateNameToNum(stateNames)))
				statecols=get_state_colors(repWhite=TRUE)[minstatenum:maxstatenum]
				statelabel_colors=get_state_colors(repWhite=TRUE)[as.numeric(mapStateNameToNum(rownames(reord.mat)))]
				pdf(outfile)
				lmat=rbind(c(4, 3, 0), c(2, 1, 0))
				lwid=c(1, 6, 1)
				lhei=c(1,6)
				cexlab=1
				par(cex.lab=cexlab)
				heatmap.3(reord.mat, lmat=lmat, lwid=lwid, lhei=lhei, main="Significant features as \n chromatin state and gene combinations", 
					density.info="none", na.rm=TRUE, na.color="white", trace="none", 
					xlab=paste(ncol(reord.mat), "significant genes"),
					labRow=rownames(reord.mat), RlabColor=statelabel_colors, labCol=rep("", length(colnames(reord.mat))),	
					cexRow=1.5, key=FALSE, dendrogram="none", col=statecols, Rowv=FALSE, Colv=FALSE, scale="none")
			}
		}, error=function(e) {
			warning(e)
			print(paste("Problem plotting ", outfile))
			save(list=c("e", "dend", "final.filled.mat", "heatcols", "reorderRows", "reorderCols", "mapStateNameToNum", "get_state_colors"), file="error.pval.rdata") 
		})
		dev.off()
		if(exists("hm") && !is.null(hm) && !is.na(hm$colInd)) {
			sigpvals.order=colnames(final.filled.mat)[hm$colInd]
			#sigpvals.dend=hm$colDendrogram
			#orig.mat=final.filled.mat
		}
	} else {
		print(too.little.data.msg); return()
	}
	write.table(sigpvals.order, row.names=FALSE, col.names=FALSE, file=paste(plotdir, filebase, ".txt", sep=""), quote=FALSE )
	sigpvals.order.file=paste(plotdir, filebase, ".rdata", sep="")
	save(list=c("sigpvals.order"), file=sigpvals.order.file)
	#save(list=c("sigpvals.order", "sigpvals.dend", "orig.mat"), file=sigpvals.order.file)
}

## ranks is rank of all p.vals
## split.mat is 2 rows, first row is gene name, 2nd row is chrstate
plot.ecdf <- function(ranks, split.mat, statevec, xmax, mainlabel, rawnum=FALSE)
{
	if(xmax>=1) {
		numstates=length(statevec)
		colors=get_state_colors()[as.character(statevec)]
		## make any white colors black so we can see it
		colors[which(colors=="#FFFFFF")] = "#000000"
		annot.info=get_annot.info()
		statenames=annot.info[as.character(statevec),"MNEMONIC"]
					
		maxes=c()
		## for each chromatin state
		for(chrstate in statevec) {
			## get features that match the chromatin state
			keep.indices=which(split.mat[2,]==chrstate)
			## get the corresponding ranks
			curr.ranks=ranks[keep.indices]
			curr.ranks=sort(curr.ranks, decreasing=FALSE)
			## get the ECDF for this chromatin state
			fn=get.state.ecdf(curr.ranks, length(split.mat[1,]), rawnum)
			maxes=append(maxes, fn(xmax))
		}
		## determine the y-axis maximum 
		ymax=max(maxes)
		for (chrstate in statevec) {
			curr.color=colors[as.character(chrstate)]	
			keep.indices=which(split.mat[2,]==chrstate)
			curr.ranks=ranks[keep.indices]
			curr.ranks=sort(curr.ranks, decreasing=FALSE)
			fn=get.state.ecdf(curr.ranks, xmax, rawnum)
			plot(fn,do.points=FALSE, col=curr.color, verticals=TRUE, ylim=c(0,ymax), xlim=c(0,xmax), main="")
			
			### TODO: ADD TICK MARKS BASED ON % GENOME OR % GENE BODY
			#ticks=blahblahblah
			#axis(side=2, at=ticks)
			#for(state in statevec) {
			#	abline(a=intercept, b=slope, col=color, lty="dashed")
			#}

			par(new=TRUE)
		}
		legend("topleft", legend=statenames, fill=colors)
		title(main=mainlabel)
	}
}



## helper function for mk.all.plots
## makes barplots of chrstate frequency and gene frequency for top 100, 1000, 10k features
make.barplots <- function(p.vals, genenames, index.order, statevec, plot.subdir, prefix) {
	statevec=as.numeric(statevec)
	numstates=length(statevec)
	chrstates=rep(statevec, length(p.vals)/numstates)
	genenames.feats=rep(genenames, each=numstates)
	sorted.pvals=p.vals[index.order]
	sorted.chrstates=chrstates[index.order]
	sorted.genenames=genenames.feats[index.order]
	#sorted.pv.file=paste("rdata/", subdir, "sorted.pvals.rdata", sep="")
	#save(sorted.pvals, file=sorted.pv.file)
	for(numtopfeats in c(100, 1000, 10000)) {
		#print(paste("Doing top ", numtopfeats, " feats", sep=""))
		pdf(paste(plot.subdir, "top", numtopfeats, prefix, "chrstates.pdf", sep=""))
		ordered.states=as.numeric(sorted.chrstates[1:numtopfeats])
		t=tabulate(ordered.states, nbins=max(statevec))
		names(t)=1:max(statevec)
		chrstate.freqs=t[statevec]
		barplot(chrstate.freqs)
		dev.off()
	
		gene.counts=table(sorted.genenames[1:numtopfeats])
		gene.freqs=tabulate(gene.counts)
		names(gene.freqs)=1:max(gene.counts)
		pdf(paste(plot.subdir, "top", numtopfeats, prefix, "genecounts.pdf", sep=""))
		barplot(gene.freqs)
		dev.off()
	}
}

draw.fig1.lines=function(horiz.ints, vert.ints, linewidth) {
	abline(h=horiz.ints, lwd=linewidth); 
	abline(v=vert.ints, lwd=linewidth)
}

## For figure 1
## saves to figure1.pdf
make.subsetfeat.heatmaps <- function(matrix.a, matrix.b, a.label, b.label, sig.feats, p.vals, base.plotdir, subdir, statevec, matched, model) {
	narm=TRUE
	type="feats"
	suffix=""
	if(matched) {
		suffix="_matched"
	}

	library("ggplot2")
	pval.outfile="fig1.pval.pdf"
	dd=with(density(p.vals), data.frame(x, y))
	ggplot(data = dd, mapping = aes(x = x, y = y)) +
		geom_line(color="black") + layer(data = dd, mapping = aes(x=ifelse(x < .05,x,.05), y=y), geom = "area", geom_params=list(fill="red",alpha=.3)) +
		scale_y_continuous(limits = c(0,max(dd$y)), name="Density") +
		scale_x_continuous(name="Corrected p-values") +
		geom_vline(aes(xintercept=.05), color="red", linetype="dashed")
	ggsave(pval.outfile)
	
	print("Plotting heatmap of some feats...")	
	## Pick some sig feature values
	sigmat.a.feats=t(matrix.a[sig.feats, ])                                     
	sigmat.b.feats=t(matrix.b[sig.feats, ])
	##sigmat.a.feats has celltypes as rows and features as columns	
	## combine as rows because we have celltypes as rows
	fullsigmat.feats=rbind(sigmat.a.feats, sigmat.b.feats)
	keep.sigfeats=sig.feats[get_order(fullsigmat.feats, 2, TRUE)][1:10]

	## find nonsig feats
	randfeats=rownames(matrix.a)[1:100]
	nonsigfeats=randfeats[which(!(randfeats %in% sig.feats))]
	nonsigmat.a.feats=t(matrix.a[nonsigfeats[c(1:20, 41:60)], ])
	nonsigmat.b.feats=t(matrix.b[nonsigfeats[c(1:20, 41:60)], ])
	fullnonsigmat.feats=rbind(nonsigmat.a.feats, nonsigmat.b.feats)

	fullmat.a=cbind(sigmat.a.feats[, keep.sigfeats], nonsigmat.a.feats)
	fullmat.b=cbind(sigmat.b.feats[, keep.sigfeats], nonsigmat.b.feats)
	fullmat=cbind(fullsigmat.feats[, keep.sigfeats], fullnonsigmat.feats)

	outfile=paste("figure1.pdf")
	#if (is.matrix(sigmat.a) && is.matrix(sigmat.b) && dim(sigmat.a)[1]>1 && dim(sigmat.b)[1]>1 && dim(sigmat.a)[2]>1 && dim(sigmat.b)[2]>1) {
		
		## pick row ordering for each group separately
		print("Ordering rows...")
		row1.ord=get_order(fullmat.a, 1, narm)
		row1.names=rownames(fullmat.a)[row1.ord]
		row2.ord=get_order(fullmat.b, 1, narm)
		row2.names=rownames(fullmat.b)[row2.ord]

		## combine and allow column reordering
		mat.to.plot=fullmat[append(row1.names, row2.names),]

		print("Generating labels...")
		## generating colors
		min=min(mat.to.plot, na.rm=TRUE)
		max=max(mat.to.plot, na.rm=TRUE)
		b1=seq(from=min, to=0, by=.01)
		b2=seq(from=0, to=max, by=.01)
		heatcol1=colorRampPalette(c("blue", "white"))(length(b1)-1)
		heatcol2=colorRampPalette(c("white", "red"))(length(b2)-1)
		heatcolors=c(heatcol1, heatcol2)
		b=c(b1, b2[2:length(b2)])

		## get chromatin states associated with feats
		## need string of chromatin states so we index by string of number, instead of direct number index, for statecolors
		allcolors=get_state_colors()
		statecolors=allcolors[as.character(statevec)]
		names(statecolors)=as.character(statevec)

		## row color labels must be matrix now
		r.col.lab=as.matrix(t(append(rep("dodgerblue4", length(row1.names)), rep("deepskyblue", length(row2.names)))))
		rownames(r.col.lab)="Grouping"

		## column color labels as chromatin state colors
		## start out with color as default black, then fill in for each state
		colstates = as.character(get_chrstates_only(colnames(mat.to.plot)))
		## apply relevant colors
		chromatin.labels=rep("black", length(colstates))
		for(currstate in unique(colstates)) {
			chromatin.labels[which(colstates==currstate)]=statecolors[currstate]
		}
		c.col.lab=matrix(chromatin.labels, ncol=1)
		
		## label rows as either group a or group b
		group.cols=c("dodgerblue4", "deepskyblue")
		## Make matrix of column labellings
		all.col.legends=c(a.label, b.label )
		all.col.fill=c(group.cols)

		## make rowname labels
		rowlabeltextcolors=get_ep_colors(rownames(mat.to.plot))
		rowlabels=get_ep_names(rownames(mat.to.plot))
		
		## make rowside labels of group names
		row.side.labels=c(a.label, b.label)
		## calc positions of group labels
		total=nrow(mat.to.plot)
		pos1=total-length(row1.names)/2
		pos2=length(row2.names)/2
		row.side.fracs=c(pos1/total, pos2/total)
				    
		## setup plotting layout
		lmat=rbind(c(4, 3, 0), c(2, 1, 0))
		lwid=c(1, 6, .1)
		lhei=c(0.75,7)

		## set up drawing of lines
		linewidth=5
		horiz.ints=c(length(row2.names)+.5)
		vert.ints=c(10.5)
		print(c(horiz.ints, linewidth, vert.ints))
		print(draw.fig1.lines)

		save(list=ls(all=TRUE), file="figure1.rdata")
		print("Plotting heatmap...")
		## actually plot matrix
		pdf(outfile, width=14, height=14)
		hm<-heatmap.3(mat.to.plot, 
			lmat=lmat, lwid=lwid, lhei=lhei,
			Colv=FALSE, Rowv=FALSE, dendrogram="none", col=heatcolors, breaks=b, 
			margins=c(2, 20), keysize=.5, scale="none", trace="none", cexRow=1,
			ylab=paste(nrow(mat.to.plot), "Epigenomes"), xlab=paste(ncol(mat.to.plot), "features"),
			RowSideColors=r.col.lab, RlabColor=rowlabeltextcolors, labRow=rowlabels, RowSideLabels=row.side.labels, RowSideFracs=row.side.fracs, cexRowSideLabels=c(2,2),
			#ColSideColors=c.col.lab, labCol=rep("", length(colnames(mat.to.plot))),
			#add.expr=draw.fig1.lines(horiz.ints, vert.ints, linewidth)
			add.expr={abline(h=horiz.ints, lwd=linewidth); abline(v=vert.ints, lwd=linewidth)}
			)
		legend("topright", legend=all.col.legends, fill=all.col.fill, border=FALSE, bty="n", y.intersp=0.7, cex=0.5)
		dev.off()
	#}
}

## matrix.a and matrix.b should have feature values in matrix with features as rows and celltypes as columns
make.sigfeat.heatmaps <- function(matrix.a, matrix.b, a.label, b.label, sig.feats, base.plotdir, subdir, matched, model, filebase) {
	if(!file.exists(paste(base.plotdir, subdir, sep=""))) {
		dir.create(paste(base.plotdir, subdir, sep=""), recursive=TRUE)
	}
	print(paste("Plotting", filebase, "feature heatmaps..."))	
	#save(matrix.a, file="tmp/mata.rdata")
	#save(matrix.b, file="tmp/matb.rdata")
	## plot based on feature values
	if(length(sig.feats) >0) {
		if(length(sig.feats)==1) {
			celltypes.a=colnames(matrix.a)
			celltypes.b=colnames(matrix.b)
			## if only one sig feat, R turns it into a vector, so we must turn it in to a matrix with the right dimensions again
			sigmat.a.feats=matrix(matrix.a[sig.feats, ], ncol=length(sig.feats), dimnames=list(celltypes.a, sig.feats))
			sigmat.b.feats=matrix(matrix.b[sig.feats, ], ncol=length(sig.feats), dimnames=list(celltypes.b, sig.feats))	
			
		} else {
			sigmat.a.feats=t(matrix.a[sig.feats, ])                                     
			sigmat.b.feats=t(matrix.b[sig.feats, ])
		}
		##sigmat.a.feats has celltypes as rows and features as columns	
		## combine as rows because we have celltypes as rows
		fullsigmat.feats=rbind(sigmat.a.feats, sigmat.b.feats)
		## plot based on rankings (within each feature)
		sigmat.a.ranks=apply(sigmat.a.feats, 2, rank)
		sigmat.b.ranks=apply(sigmat.b.feats, 2, rank)
		fullsigmat.ranks=apply(fullsigmat.feats, 2, rank)
	
		save(list=ls(all=TRUE), file="tmp.rdata")
		for(run in c(1)) {
		#for(run in 1:2) {
			if (run==1) {
				sigmat.a=sigmat.a.feats
				sigmat.b=sigmat.b.feats
				fullsigmat=fullsigmat.feats
				type="feats"
			} else {
				sigmat.a=sigmat.a.ranks
				sigmat.b=sigmat.b.ranks
				fullsigmat=fullsigmat.ranks
				type="ranks"
			}

			suffix=""
			if(matched) {
				suffix="_matched"
			}
			outfile=paste(base.plotdir, subdir, filebase, "_", type, "_plot", suffix, ".pdf", sep="")
			featorder.file=paste(base.plotdir, subdir,filebase, "_", type,"_featorder", suffix, ".txt", sep="")
			celltypeorder.file=paste(base.plotdir, subdir, filebase,"_", type,"_celltypeorder", suffix, ".txt", sep="")
			
			dendfile=paste(base.plotdir, subdir, filebase,"_", type,"_dend", suffix, ".rdata", sep="")
			d1file=paste(base.plotdir, subdir, filebase,"_", type,"_d1", suffix, ".rdata", sep="")
			d3file=paste(base.plotdir, subdir, filebase,"_", type,"_d3", suffix, ".rdata", sep="")
			mat.file=paste(base.plotdir, subdir, filebase,"_", type,"_mat.to.plot", suffix, ".rdata", sep="")
			groupings.file=paste(base.plotdir, subdir, filebase,"_", type,"_colname_groupings", suffix, ".rdata", sep="")
			
			if(length(sig.feats)==1){
				celltype.annots=append(paste(celltypes.a, "_", a.label, sep=""), paste(celltypes.b, "_", b.label, sep=""))
				plot_vector_as_hm(outfile, fullsigmat,  celltype.annots, colnames(fullsigmat))
				write.table(sig.feats, file=featorder.file, quote=FALSE, row.names=FALSE)
				write.table(rownames(fullsigmat), file=celltypeorder.file, quote=FALSE, row.names=FALSE)
				mat.to.plot=fullsigmat
				save(mat.to.plot, file=mat.file)
				next()	
			}
			
			if (is.matrix(sigmat.a) && is.matrix(sigmat.b) && dim(sigmat.a)[1]>1 && dim(sigmat.b)[1]>1 && dim(sigmat.a)[2]>1 && dim(sigmat.b)[2]>1) {
				## first pick column ordering from all celltypes
				print("Separating into three groups and ordering columns and rows...")
				narm=TRUE
	
				## pick row ordering for each group separately
				row1.ord=get_order(sigmat.a, 1, narm)
				row1.names=rownames(sigmat.a)[row1.ord]
				
				row2.ord=get_order(sigmat.b, 1, narm)
				row2.names=rownames(sigmat.b)[row2.ord]
	
				## pick column ordering overall
				colorder=get_order(fullsigmat, 2, narm)
				## combine both groups of rows and allow column reordering within three column groups
				ordered.mat=fullsigmat[append(row1.names, row2.names), colorder]
				
				## find if genes that show up more than once	
				genes<-get_genes_only(colnames(ordered.mat))
				names(genes)=colnames(ordered.mat)
				gene_table <- table(genes)
			 	repeat_genes <- names(gene_table)[which(gene_table>=2)]
			 	norepeats=FALSE
				if(length(repeat_genes)==0) {norepeats=TRUE}	
				save(list=ls(all=TRUE), file="tmp2.rdata")	
				## without matching or with no repeat genes, just directly order columns of matrix
				if(!matched || norepeats) {
					if(ncol(ordered.mat)==1) { dend.featnames=colnames(ordered.mat) } else {
						dend=load_or_get_dend(dendfile, ordered.mat, 2, narm, "dend")
						dend.featnames=colnames(ordered.mat)[order.dendrogram(dend)]
					}
					mat.to.plot=ordered.mat[, dend.featnames]
				} else {
					## there are repeat genes, so split into 3 groups
			 		## take their first instance (based on genes in col.names) and put them in group 1
			 		first_inds_unordered=match(repeat_genes, genes)
			 		group1_feats=names(genes)[sort(first_inds_unordered, decreasing=FALSE)]
			 		group1_genenames_inorder=get_genes_only(group1_feats)

			 		## find their matches for group 2, order same 
			 		all_feats=colnames(ordered.mat)
			 		leftover_feats=all_feats[!(all_feats %in% group1_feats)]
			 		names(leftover_feats)=get_genes_only(leftover_feats)
			 		## find other indices for each gene (in order)
			 		## unlist results for ones repeated more than twice, just keep them in same order
			 		group2_feats_ordered_by_genenames=unlist(sapply(group1_genenames_inorder, function(x) { leftover_feats[which(names(leftover_feats)==x)]} ), use.names=FALSE)

			 		unique_genes=names(gene_table)[which(gene_table==1)]
			 		unique_inds=match(unique_genes, genes)
					group3_feats=names(genes)[sort(unique_inds, decreasing=FALSE)]
					
					print("Ordering columns within three groups...")
				    group1_mat=ordered.mat[,group1_feats]
                    if(length(group1_feats)==1) {
                        ## put in this clause because it automatically converts one column into a vector
                        group1_mat=matrix(ordered.mat[,group1_feats], dimnames=list(rownames(ordered.mat), group1_feats))
                    }
                    if(ncol(group1_mat)==1) { final.group1.featnames=colnames(group1_mat) } else {
                    	d1=load_or_get_dend(d1file, group1_mat, 2, narm, "d1")
						#hclust=hclust(dist(t(group1_mat), method="euclidean"), method="complete")
						#dend=as.dendrogram(hclust)
						final.group1.featnames=colnames(group1_mat)[order.dendrogram(d1)]
					}
					final.group1.genenames=get_genes_only(final.group1.featnames)

					## apply groupings and orderings to group2 based on group 1
					cnames=colnames(ordered.mat)
					final.leftover.feats=cnames[!(cnames %in% final.group1.featnames)]
					final.leftover.genenames=get_genes_only(final.leftover.feats)
					final.group2.featnames=unlist(sapply(final.group1.genenames, function(x) {final.leftover.feats[which(final.leftover.genenames==x)]}))

					## now get dendrogram and new ordering for subgroup 3
					group3_mat=ordered.mat[,group3_feats]
					if(length(group3_feats)==1) {
						## put in this clause because it automatically converts one column into a vecotr
						group3_mat=matrix(ordered.mat[,group3_feats], dimnames=list(rownames(ordered.mat), group3_feats))
					}
					if(ncol(group3_mat)==0) { 
						final.group3.featnames=c() 
					} else if(ncol(group3_mat)==1) { 
						final.group3.featnames=colnames(group3_mat) 
					} else {
						#save(list=ls(all=TRUE), file="tmp.rdata")
						d3=load_or_get_dend(d3file, group3_mat, 2, narm, "d3")
						#hclust=hclust(dist(t(group3_mat), method="euclidean"), method="complete")
						#dend=as.dendrogram(hclust)
						final.group3.featnames=colnames(group3_mat)[order.dendrogram(d3)]
					}
					## make labelings of groups
					genegroup.legends=c("Group 1", "Group 1 matched", "Unique")
					genegroup.colors=c("mistyrose1", "mistyrose3", "mistyrose4")
					genegroup.labels=c(rep(genegroup.colors[1], length(final.group1.featnames)), rep(genegroup.colors[2], length(final.group2.featnames)), rep(genegroup.colors[3], length(final.group3.featnames)))
					
					ordered.colnames=c(final.group1.featnames, final.group2.featnames, final.group3.featnames)
					print(dim(ordered.mat))
					print(length(ordered.colnames))
					matched_ordered.mat=ordered.mat[,ordered.colnames]
					groupings=c(rep(1, length(final.group1.featnames)), rep(2, length(final.group2.featnames)), rep(3, length(final.group3.featnames)))
					names(groupings)=colnames(matched_ordered.mat)
					save(groupings, file=groupings.file)
					mat.to.plot=matched_ordered.mat
				}
				
				print("Generating final heatmap...")
				## row labels are grouping (group 1 vs group 2) 
				epcolor.labels=get_ep_colors(rownames(mat.to.plot))
				group.labels=append(rep("dodgerblue4", length(row1.names)), rep("deepskyblue", length(row2.names)))
				rlab=rbind(epcolor.labels, group.labels)
				rownames(rlab)=c("Coloring", "Grouping")

				heatcolors=colorRampPalette(c("blue", "white", "red"))(100)
		
				#metric=strsplit(subdir, "/")[[1]][1]
				allcolors=get_state_colors()
				statevec=unique(get_chrstates_only(colnames(mat.to.plot)))
				statecolors=allcolors[as.character(statevec)]
				names(statecolors)=statevec
	
				print("Generating chromatin state labeling...")
				## get column labels of chromatin state, background vals, and pvals
				## get chromatin states associated with feats
				## need string of chromatin states so we index by string of number, instead of direct number index, for statecolors
				colstates = as.character(get_chrstates_only(colnames(mat.to.plot)))
				## apply relevant colors
				chromatin.labels=rep("black", length(colstates))
				for(currstate in unique(colstates)) {
					chromatin.labels[which(colstates==currstate)]=statecolors[currstate]
				}

				print("Generating labellings...")
				## make rowname labels
				rowlabeltextcolors=get_ep_colors(rownames(mat.to.plot))
				rowlabels=get_ep_names(rownames(mat.to.plot))
				
				group.cols=c("dodgerblue4", "deepskyblue")
				if(!matched || norepeats) {
					## Make matrix of column labellings
					clab=as.matrix(chromatin.labels)
					colnames(clab)=c("Chromatin_state")
					all.col.legends=c(a.label, b.label)
					all.col.fill=c(group.cols)
				} else {
					if(length(genegroup.labels)!=length(chromatin.labels)) {error("Length of gene labels and chromatin labels don't match.") }
					clab=cbind(chromatin.labels, genegroup.labels)
					colnames(clab)=c("Chromatin_state", "Gene type")
					all.col.legends=c(a.label, b.label, "",genegroup.legends)
					all.col.fill=c(group.cols, "white", genegroup.colors)
				}

				print("Plotting matrix...")
				save(list=ls(all=TRUE), file="fig1.rdata")
				## actually plot matrix
				pdf(outfile, width=14, height=14)
				hm<-heatmap.3(mat.to.plot, RlabColor=rowlabeltextcolors, margins=c(2, 20), Colv=FALSE, Rowv=FALSE, dendrogram="none", keysize=.75, ylab=paste(nrow(mat.to.plot), "Epigenomes"), xlab=paste(ncol(mat.to.plot), "features"),scale="none", trace="none",RowSideColors=rlab, labCol=rep("", length(colnames(mat.to.plot))), col=heatcolors, ColSideColors=clab, labRow=rowlabels)
				legend("topright", legend=all.col.legends, fill=all.col.fill, border=FALSE, bty="n", y.intersp=0.7, cex=0.5)
				dev.off()

				## save feature order of columns
				col.permutation=hm$colInd
				feat.cols=colnames(mat.to.plot)[col.permutation]
				write.table(feat.cols, file=featorder.file, quote=FALSE, row.names=FALSE)
				write.table(rownames(mat.to.plot), file=celltypeorder.file, quote=FALSE, row.names=FALSE)
				save(mat.to.plot, file=mat.file)
			} else {
				print(paste("Not enough significant figures for ", outfile, sep=""))
			}
		} 
		rm(sigmat.a.feats)
		rm(sigmat.b.feats)
		rm(fullsigmat.feats)
		rm(sigmat.a.ranks)
		rm(sigmat.b.ranks)
		rm(fullsigmat.ranks)
	}
}


mk.all.pval.plots <- function(sig.pvals, base.plotdir, subdir, model) {
	## turn off for debugging purposes
	heatmaps=TRUE
	ecdf=TRUE
	plot.subdir=paste(base.plotdir, subdir, sep="")
	dir.create(plot.subdir, recursive=TRUE, showWarnings=FALSE)
	num.sig = length(sig.pvals)	
	
	## replace NA's in pval matrix with 1
	indices=which(is.na(sig.pvals))
	sig.pvals[indices]=1
	
	statevec=unique(get_chrstates_only(names(sig.pvals)))
	genenames=unique(get_genes_only(names(sig.pvals)))
	numstates=length(statevec)

	sigpval.mat=matrix(rep(NA, numstates*length(genenames)), nrow=numstates, dimnames=list(statevec, genenames))
	for(sigfeat in names(sig.pvals)) {
		currstate=get_chrstates_only(sigfeat)
		currgene=get_genes_only(sigfeat)
		sigpval.mat[currstate, currgene] = sig.pvals[sigfeat]
	}

	## heatmap of matrix of pvals has chr states as rows and genes as columns
	print("Generating matrix of p-values...")
	## identify orderings of features
	index.rank=get.ranks.with.ties.for.nas(as.vector(sigpval.mat))

	if(heatmaps) {	
		## make color scale blue to white for significant pvals, then white for the rest
		## set any non-sig pvals to 1
		colfunc=colorRampPalette(c("blue", "white"))
		heatcols=colfunc(100)
		## plot significant only 
		print("Plotting heatmap of significant feature pvals...")
		plot.top.n(sigpval.mat, index.rank, num.sig, TRUE, TRUE, plot.subdir, "sig.pvals", heatcols, model, TRUE)
	}
	
	################ EMPRICIAL CUMULATIVE DISTRIBUTION FUNCTION PLOTS #######################
	if(ecdf) {
		## for chrstates 
		print("Plotting empirical cumulative distribution plot...")
		ranks=rank(sig.pvals)
		split.mat=get_split_mat(names(sig.pvals))


		## plot sig only
		ecdf.file=paste(plot.subdir, "sig.ecdf.chrstates.pdf", sep="")
		pdf(ecdf.file)
		mainlabel="State ECDF for sig pvals"
		plot.ecdf(ranks, split.mat, statevec, num.sig, mainlabel, rawnum=TRUE)
		dev.off()
	}
}
