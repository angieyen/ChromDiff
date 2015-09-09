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
suppressMessages(library(methods, quietly=TRUE))
suppressMessages(library(utils, quietly=TRUE))

## if no data, return NA 
my.which.max <- function(x) {
	if(length(which.max(x))==0) {
		return(NA)
	} else {
		return(which.max(x))
	}
}
my.which.min <- function(x) {
	if(length(which.min(x))==0) {
		return(NA)
	} else {
		return(which.min(x))
	}
}

#get_corrected_logexp=function(group1, group2, genes) {
#	orig.data=get_logexp(c(group1, group2), genes)
#	group1.mean=mean(orig.data[group1, genes], na.rm=TRUE)
#	group2.mean=mean(orig.data[group2, genes], na.rm=TRUE)
#	final.mat=rbind(orig.data[group1, genes]-group1.mean, orig.data[group2, genes] - group2.mean)
#	return(final.mat)
#}
get_random_logexp=function(celltype.order, numgenes) {
	logexp.mat=get_all_logexp()
	random.cols=sample(1:ncol(logexp.mat), numgenes)
	random.genes=colnames(logexp.mat)[random.cols]
	return(get_logexp(celltype.order, random.genes))
}
get_other_genes=function(genes.to.skip) {
	logexp.mat=get_all_logexp()
	allgenes=colnames(logexp.mat)
	genes.to.keep=allgenes[which(!(allgenes %in% genes.to.skip))]
	return(genes.to.keep)
}
get_other_logexp=function(celltype.order, genes.to.skip) {
	genes.to.keep=get_other_genes(genes.to.skip)
	return(get_logexp(celltype.order, genes.to.keep))
}

## get log expression for requested ensembl genes and celltypes
get_logexp=function(celltype.order, ensembl.geneorder) {
    logexp.mat=get_all_logexp()
    logexp.mat=fill_nas(logexp.mat, celltype.order, ensembl.geneorder)
    return(logexp.mat)
}
fill_nas = function(mat, out.rownames, out.colnames) {
    default=NA          
    ## fill in colnames
    na.colnames=out.colnames[which(!out.colnames %in% colnames(mat))]
	fill.cols=matrix(rep(default, length(na.colnames)*nrow(mat)), nrow=nrow(mat), dimnames=list(rownames(mat), na.colnames))
	mat=cbind(mat, fill.cols)
                	
	na.rownames=out.rownames[which(!out.rownames %in% rownames(mat))]
	fill.rows=matrix(rep(default, length(na.rownames)*ncol(mat)), ncol=ncol(mat), dimnames=list(na.rownames, colnames(mat)))
	mat=rbind(mat, fill.rows)
	result.mat=mat[out.rownames, out.colnames]
	return(result.mat)
}


## given a matrix, generate the dendrogram for it for dimension "dim", save it in savefile 
load_or_get_dend <- function(savefile, mat, dim, narm, dend.str) {
    if(any(dim(mat)==0)) { print("Warning: matrix passed in has dimension of 0"); return();}
	if(dim(mat)[dim]==1) { stop(paste("Can not get dend for matrix with length of 1 for dimension", dim)) }
	
    ## if can load everything from saved file, do so and return
    success=FALSE
    if(file.exists(savefile)) {
        print(paste("Loading", savefile, "..."))
        load(savefile)
        if(exists(dend.str)) {
            success=TRUE
        }
		## if file exists but has the wrong data stored, delete the file
        if(!success) {
        	file.remove(savefile)
        }
    }
   

    # if could not load information, re-calculate now
    if(!success) {
    	#save(list=ls(all=TRUE), file="tmp.rdata")
    	mat=as.matrix(mat)
        print("Recalculating dend...")
        dend=get_dend(mat, dim, narm)
		assign(dend.str, dend)
        print("Saving dend info...")
		parentdir=dirname(savefile)
		dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
		save(list=c(dend.str), file=savefile)
    }
    return(get(dend.str))
}

comp.dens.plot = function(savefile, group1.vals, group2.vals, group1, group2) {
    suppressMessages(library("ggplot2"))
    allvals=c(group1.vals, group2.vals)
    labels=c(rep(group1, length(group1.vals)), rep(group2, length(group2.vals)))
    df=data.frame(value=allvals, label=labels)
    p=ggplot(df) + geom_density(aes(x = value, group=label, colour = label)) + ggtitle("Densities from a kernel density estimator")
    ggsave(savefile)
}
rm.nas = function(vec) {
    return(vec[which(!is.na(vec))])
}

## get appropriate dend via hclust, ordered by row or column mean (for dim 1 or 2, respectively)
get_dend <- function(mat, dim, narm) {
	if(dim==1) {
		transformfunc<-function(x) {return(x)}
		meansfunc=rowMeans
	} else if (dim==2) {
		transformfunc<-function(x) {return(t(x))}
		meansfunc=colMeans
	} else {
		stop(paste("unrecognized dimension ", dim, sep=""))
	}
	success=tryCatch({
		reorder(as.dendrogram(hclust(dist(transformfunc(mat)))), meansfunc(mat, na.rm=narm))
		TRUE
	}, error=function(e) {
		FALSE
	})

	if(!success) {
		## replace NAs with 0's
		mat[which(is.na(mat))]=0
	}
	dend=reorder(as.dendrogram(hclust(dist(transformfunc(mat)))), meansfunc(mat, na.rm=narm))
	return(dend)
}

get_order <- function(mat, dim, narm) {
	return(order.dendrogram(get_dend(mat, dim, narm)))
}


gen_colors <- function(min, med, max) {
	range=max-min
	scalingfactor=100/range
	nbw=(floor((med-min)*scalingfactor))
	nwr=floor((max-med)*scalingfactor)
	bwcols=colorRampPalette(c("blue", "white"))(nbw)
	wrcols=colorRampPalette(c("white", "red"))(nwr)
	return(append(bwcols, wrcols))
}

get_features_from_genes = function(genes, featurenames) {
	corresponding_genes=get_genes_only(featurenames)
	assoc.feats=names(corresponding_genes[which(corresponding_genes %in% genes)])
	return(assoc.feats)
}
pick_recurring_gene_feats=function(considered.feats, cutoff) {
	genenames_count=table(get_genes_only(names(considered.feats)))
	ord_genenames_count=sort(genenames_count, decreasing=TRUE, method="shell")
	runsum=0
	genes_to_use=c()
	for(ind in 1:length(ord_genenames_count)) {
		if(runsum>cutoff) {
			break
		}
		runsum=runsum+ord_genenames_count[ind]
		genes_to_use=append(genes_to_use, names(ord_genenames_count)[ind])
	}
	## get features associated with those genes
	feats_to_keep=get_features_from_genes(genes_to_use, names(considered.feats))[1:cutoff]
	return(feats_to_keep)
}

## returns matrix with features as rows and cell types as columns
get_residuals <- function(model, metric, dependentvariable, no_covariate_correction) { 
	if(dependentvariable=="special") {
		dependentvariable="age"
	}else if (dependentvariable=="specialgi") {
		dependentvariable="anatomy"
	}
	## go through each celltype, read in corrected values (calculate corrected values if they don't exist yet)
	matrix.vec=c()
	prefix=get_prefix(metric, model)
	metricdir=get_metric_subdir(metric)
	suffix=get_suffix(metric)
	
	## check if dependent variable is one of the covariates we correct for
	if(toupper(dependentvariable) %in% get_covariates_to_correct()) {
		varsuffix=paste(".skip", dependentvariable, sep="")
	} else {
		varsuffix=".skipnone"
	}

	residmat.rdatafile=paste("rdata/", model, "/", metricdir, "residmat", varsuffix, suffix, ".rdata", sep="")

	success=FALSE
	if(file.exists(residmat.rdatafile) && file.info(residmat.rdatafile)$size>0) {
		print(paste("Loading residuals from ", residmat.rdatafile, "...", sep=""))
		tryCatch( {load(residmat.rdatafile); 
					if(nrow(resids)<ncol(resids)) {resids=t(resids)};
					if(all(is.na(resids)==FALSE)) {success=TRUE}
					}, 
			error=function(e) {print("Error reading residual matrix file...")}) 
	
	}
	#print(colnames(resids))
	if(!success) {
		featnames=get_featnames(model, metric)		
		allcelltypes.rdatafile=paste("rdata/", model, "/", metricdir,  "all" , suffix, ".rdata", sep="")
		matLoaded=FALSE
		if(file.exists(allcelltypes.rdatafile) && file.info(allcelltypes.rdatafile)$size>0) {
			print(paste("Loading feature values for all celltypes from", allcelltypes.rdatafile, "..."))
			tryCatch( {load(allcelltypes.rdatafile); 
				## make sure that necessary variables exists, and make sure there were no errors with the matrix (there should be no NA elements)
				if(exists("mat") && (metric!="perc" || exists("totalcounts")) && (all(is.na(mat)==FALSE))) {matLoaded=TRUE}},
				error=function(e){print("Error reading all epigenomes files...")})
		}
		if(!matLoaded) {
			firstcounts=TRUE
			celltypes.list=get_all_celltypes()
			print("Loading celltype data from each file...")
			validcelltypes=c()
			for (ind in 1:length(celltypes.list)) {
				vecLoaded=FALSE
				sharedrdatafile=paste("rdata/", model, "/", metricdir,  celltypes.list[ind], suffix, ".rdata", sep="")
				datafile=paste(prefix, celltypes.list[ind], suffix, ".txt", sep="")
				## load or read in data of feature values for each celltype
				if (file.exists(sharedrdatafile) && file.info(sharedrdatafile)$size>0) {
					print(paste0("Loading rdata file for ", celltypes.list[ind], "..."))
					load(sharedrdatafile)
					if(all(is.na(curr.vec)==FALSE)) {
						vecLoaded=TRUE
					}
				} 
				if (!vecLoaded) {
					if (file.exists(datafile) && file.info(datafile)$size>0) {
						print(paste0("Reading data file for ", celltypes.list[ind], "..."))
						curr.vec=scan(datafile, quiet=TRUE)
						parentdir=dirname(sharedrdatafile)
						dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
						save(curr.vec, file=sharedrdatafile)
					} else { ## feature values txt file did not exist for this epigenome
						warning(paste0(datafile, " not found. Skipping ", celltypes.list[ind], "..."))
						next
					}
				}
				validcelltypes=append(validcelltypes, celltypes.list[ind])
				if(metric=="perc") {
					countsfile=paste(prefix, celltypes.list[ind], "_counts.txt", sep="")
					newtotalcounts=scan(countsfile, quiet=TRUE)
					if(firstcounts) {
						totalcounts=newtotalcounts
						firstcounts=FALSE
					} else {
						if(!all(totalcounts==newtotalcounts)) {
							stop(paste("Total counts of", celltypes.list[ind], "does not agree with previous counts"))
						}
					}
				}
				if(length(curr.vec)!=length(featnames)) {
					stop(paste("Wrong number of feature numbers for ", celltypes.list[ind], ": expected ", length(featnames), " but saw ", length(curr.vec), sep=""))
				}
				matrix.vec=append(matrix.vec, curr.vec)
			}
			
			if(length(matrix.vec)==0) {
				print("No state calls found for any cell type...")	
				quit()
			}
			## matrix has features as rows and cell types as columns
			mat=matrix(matrix.vec, ncol=length(validcelltypes), dimnames=list(featnames, validcelltypes), byrow=FALSE)
			parentdir=dirname(allcelltypes.rdatafile)
			dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
			validcelltypes.rdatafile=get_validcelltypesfile(model, metric, dependentvariable)
			save(validcelltypes, file=validcelltypes.rdatafile)
			if(metric=="perc") {
				save(mat, totalcounts, file=allcelltypes.rdatafile)
			} else {
				save(mat, file=allcelltypes.rdatafile)
			}
		}
		
		if(no_covariate_correction) {
			return(mat)
		}

		## covariate mat has variables as columns and cell types as rows
		##correct for all covariates EXCEPT current dependent variable
		fixed.cov.mat <- get_covariate_mat(dependentvariable)
		## pick out valid celltypes (listed in matrix of data) only
		fixed.cov.mat = fixed.cov.mat[colnames(mat), ]
		
		## invert matrix for helper function
		inputmat=t(mat)
		if(metric=="perc") {
			resids=get_residuals_from_mats(fixed.cov.mat, inputmat, totalcounts, method="log")
		} else {
			resids=get_residuals_from_mats(fixed.cov.mat, inputmat, method="lin")
		}
		## we need resids to have same dimensions as original mat, so invert again
		resids=t(resids)

		print("Saving residuals...")
		parentdir=dirname(residmat.rdatafile)
		dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
		save(resids, file=residmat.rdatafile) 
	}
	return(resids)
}

## expects the same number of rows in cov.mat and values.mat, as they are the data being corrected (epigenomes)
## cov columns factors/covariates
## values.mat columns the values (genes) being corrected
get_residuals_from_mats = function(cov.mat, values.mat, method=c("log", "lin"), totalcounts=c())  {
	if(method=="log") {
		if(length(totalcounts)!=ncol(values.mat)) {
			stop("Counts needed for logistic regression: length of totalcounts must match number of rows in values.mat")
		}
	}
	if(nrow(cov.mat)!=nrow(values.mat)) {
		stop("cov.mat and values.mat inputted into get_residuals_from_mats must have matching number of rows (epigenomes)")
	}

	## covariate matrix has celltypes as rows and factors/covariates as columns
	if(ncol(cov.mat)>=1) {
		if(method=="log") {
			formstr=" success_fail ~ cov.mat[,1] "
		} else if (method=="lin") {
			formstr=" currvals ~ cov.mat[,1] "
		}

		for(colnum in 2:ncol(cov.mat)) {
			formstr=paste(formstr, " + cov.mat[,", colnum, "]", sep="")
		}
		print(formstr)
		print("Calculating residuals...")
		
		## if any feature has any NA's, make the residual NAs
		resids=matrix(rep(0, times=nrow(values.mat)*ncol(values.mat)), ncol=ncol(values.mat), dimnames=list(rownames(values.mat), colnames(values.mat)))
		## correct each gene/region one at a time (each column)
		for(i in 1:ncol(values.mat)) {
			currvals=values.mat[,i]	
			if(length(which(is.na(currvals)))>0) { 
				result=(rep(NA, length(currvals)))
			} else {
				if(method=="log") {
					counts=rep(totalcounts[i], times=length(currvals))
					success_fail=cbind(round(currvals*counts), round((1-currvals)*counts))
					result=(resid(glm(as.formula(formstr), family="binomial"), type="deviance"))
				} else if(method=="lin"){
					result=(resid(lm(as.formula(formstr))))
				} else {
					stop('Unrecognized method for get_residuals_from_mats: should be lin or log')
				}
			}
			resids[,i]=result	
			if((i%%10000)==0) {
				print(paste("Processing feature ", i, "...", sep=""))
			}
		}
	} else {
		resids=values.mat
	}
	return(resids)
}


load.background.vals <- function(metric, featnames, model) {
	bgdir=paste("backgrounds/", model, "/", sep="")
	if(metric=="windows" || metric=="deltawindows") {
		fn="window_background"
	} else if(metric=="deltas" || metric=="perc"){
		fn="perc_background"
	} else {
		stop(paste("Unrecognized metric: ", metric))
	}
	rdatafile=paste(bgdir, fn, ".rdata", sep="")
	if(file.exists(rdatafile) && file.info(rdatafile)$size>0) {
		load(rdatafile)
	} else {
		fullfile=paste(bgdir, fn, ".txt",  sep="")
		bg.vals <- scan(fullfile)
		names(bg.vals)=featnames
		parentdir=dirname(rdatafile)
		dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
		save(bg.vals, file=rdatafile)
	}
	return(bg.vals)
}

##  http://adv-r.had.co.nz/memory.html#garbarge-collection
mem <- function() {
	bit <- 8L * .Machine$sizeof.pointer
    if (!(bit == 32L || bit == 64L)) {
        stop("Unknown architecture", call. = FALSE)
	}

	node_size <- if (bit == 32L) 28L else 56L

	usage <- gc()
    sum(usage[, 1] * c(node_size, 8)) / (1024 ^ 2)
}


comp.groups <- function(property, group.a, group.b, a.label, b.label, metrics="perc", tests="wilcox", test.corrections="fdr", model="core", no_covariate_correction=FALSE) {
	pair.label=paste(a.label, ".", b.label, sep="")
	
	## read gene names
	genenames_file=paste0("genes/", generegions_label, "/genenames.txt")
	genenames=scan(genenames_file, what=character())
	## feature information

	check_input(metrics, tests, test.corrections)
	## start doing requested test and metric combinations
	for (test in tests) {
		test.func=get_test_func(test)

		## repeat analysis for both percentage features and delta (diff from background) features
		for(metric in metrics) {

			prefix=get_prefix(metric, model)
			suffix=get_suffix(metric)
			
			full.rdatadir=get_full_rdatadir(model, test, metric)
			dir.create(full.rdatadir, recursive=TRUE, showWarnings=FALSE)

			for(correction in test.corrections) {
				pv.file=paste(full.rdatadir, pair.label, ".pvals.rdata", sep="")
				tmp.pv.file=paste(full.rdatadir, "fdr.", pair.label, ".pvals.rdata", sep="")
				corr.func=get_corr_func(correction)
					
				### Get corrected residuals (instead of raw data)
				## residuals have celltypes as columns and features as rows
				print("Getting residuals...")
				resid.mat=get_residuals(model, metric, property, no_covariate_correction)
				## pick out celltypes you are interested in
				valid_celltypes_a=group.a[which(group.a %in% colnames(resid.mat))]
				valid_celltypes_b=group.b[which(group.b %in% colnames(resid.mat))]
				matrix.a=resid.mat[,valid_celltypes_a]
				matrix.b=resid.mat[,valid_celltypes_b]

				if(dim(matrix.a)[1]!=dim(matrix.b)[1]) {
					stop("ERROR: number of features (rows) in matrices should match")
				}

				
				featnames=get_featnames(model, metric)
			
				save(list=ls(all=TRUE), file="tmp.rdata")
				pvals.success=FALSE
				## read pval from rdata file pv.file (NOTE: doesn't need correction, because they are raw p-values, but for legacy reasons, we try tmp.pv.file which uses the pv.file with "fdr" in the filename)
				if(file.exists(tmp.pv.file) && file.info(tmp.pv.file)$size>0) {
					print(paste("Loading pvals from", tmp.pv.file, "..."))
                    tryCatch( {load(tmp.pv.file); 
                        if(length(p.vals)==length(featnames)) {
                            pvals.success=TRUE
                        }}, 
                        error=function(e) {print("Error reading pval rdata file...")})
				}
				if(!pvals.success && file.exists(pv.file) && file.info(pv.file)$size>0){
					print(paste("Loading pvals from", pv.file, "..."))
					tryCatch( {load(pv.file); 
						if(length(p.vals)==length(featnames)) {
							pvals.success=TRUE
						}}, 
						error=function(e) {print("Error reading pval rdata file...")})
				}
				if(!pvals.success) {
					print("Calculating pvals...")
					## do test for each feature, remember the indices of features with significant values
					p.vals=test.func(matrix.a, matrix.b)		
					parentdir=dirname(pv.file)
					dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
					save(p.vals,file=pv.file)
					pvals.success=TRUE
				}	
				if(length(names(p.vals))==0) {
					names(p.vals)=featnames
					parentdir=dirname(pv.file)
					dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
					save(p.vals,file=pv.file)
				}
			
				
				## multiple hypothesis correction
				corrected.pvals=corr.func(p.vals)
				
				allpval.file=get_pval_file(model, metric, test, correction, property, a.label, b.label)
				print(paste("Writing all pvals to ", allpval.file, "...", sep=""))
				parentdir=dirname(allpval.file)
				dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
	 			write.table(corrected.pvals, quote=FALSE, col.names=FALSE, file=allpval.file)

				## use significant corrected p-vals< .05 OR top 2000 if more than 2000 are significant
				pvalcutoff=.05
				cutoff=10000
				all.sig.feats=names(which(corrected.pvals<pvalcutoff))
				if(length(which(corrected.pvals<pvalcutoff))>cutoff) {
					considered.feats=corrected.pvals[which(corrected.pvals<pvalcutoff)]
					sig.feats=pick_recurring_gene_feats(considered.feats, cutoff)		
				} else {
					sig.feats=names(which(corrected.pvals<=pvalcutoff))
				}
				
				## save significant features
				if(length(sig.feats)>0) {
					print(paste(length(sig.feats), "significant features found. Now visualizing..."))
					sigpval.file=get_sigpval_file(model, metric, test, correction, property, a.label, b.label)
					print(paste("Writing pvals to ", sigpval.file, "...", sep=""))
					parentdir=dirname(sigpval.file)
					dir.create(parentdir, recursive=TRUE, showWarnings=FALSE)
					write.table(corrected.pvals[sig.feats], quote=FALSE, col.names=FALSE, file=sigpval.file)
				} else {
					print("No significant features found.")
				}
			}
		}
	}
}

reorder_get_vert_ints_height=function(colnames, dend.to.reorder, heightcutoff) {
	tmp.colnames=colnames[order.dendrogram(dend.to.reorder)]
	return(get_vert_ints_height(tmp.colnames, dend.to.reorder, heightcutoff) )
}

## returns a list of clusters, with each showing first and last element of the cluster
## use height cutoff
get_vert_ints_height = function(colnames, dend, heightcutoff) {
	if(attributes(dend)$members==1) {
		vert.ints=list(c(1,1))
	} else {
		## cut dendrogram at appropriate height
		robust=cut(dend, heightcutoff)
		nclust=length(robust$lower)
		vert.ints=list()
		for(ind in 1:nclust) {
			curr.dend=robust$lower[[ind]]
			## match dendrogram members (labels) to existing order of columns in matrix (colnames)
			indices=match(labels(curr.dend), colnames)
			## if there wer any problems, print info and stop
			if(is.na(indices[1]) || !all(sort(indices)==(min(indices):max(indices)))) { 
				print(indices)
				print(curr.dend)
				print(colnames)
				save(list=ls(all=TRUE), file="error.vert.rdata")
				stop("Problem with getting vert indices: dend and colnames are not consistent");
			}
			## indices should be all consecutive
			if(!all(sort(indices)==(min(indices):max(indices)))) {stop("Current ordering of column names conflicts with dendrogram used.")}
			vert.ints[[ind]]=c(min(indices), max(indices))
			#curr.sum=curr.sum+curr.clust.size
			#vert.ints=append(vert.ints, curr.sum+.5)
		}
	}	
	return(vert.ints)
}
#get_vert_ints = function(colnames, pvclust, cutoff) {
#	print(pvclust)
#	if(nrow(pvclust$edges)==1) {
#		vert.ints=list(c(1,1))
#	} else {
#		robust=pvpick(pvclust, alpha=cutoff)
#		nclust=length(robust$clusters)
#		vert.ints=list()
##		for(ind in 1:nclust) {
#			curr.clust=robust$clusters[[ind]]
#			indices=match(curr.clust, colnames)
#			## if there wer any problems, print info and stop
#			if(is.na(indices[1]) || !all(sort(indices)==(min(indices):max(indices)))) { 
#				print(indices)
#				print(curr.clust)
#				print(colnames)
#				save(list=ls(all=TRUE), file="error.vert.rdata")
#				stop("Problem with getting vert indices: pvclust and colnames are not consistent");
#			}
#			## indices should be all consecutive
#			if(!all(sort(indices)==(min(indices):max(indices)))) {stop("Current ordering of column names conflicts with dendrogram used.")}
#			vert.ints[[ind]]=c(min(indices), max(indices))
#			#curr.sum=curr.sum+curr.clust.size
#			#vert.ints=append(vert.ints, curr.sum+.5)
#		}
#	}	
#	return(vert.ints)
#}

## get index of last instance of each element in tobematched in tosearch
match_last <- function(tobematched, tosearch) {
	return(unlist(sapply(tobematched, function(tomatch){return(max(which(tomatch==tosearch)))})))

}
order.list.by.first=function(mylist) {
    first.elts=sapply(mylist, function(x) x[1])
	reord=order(first.elts)
	new.list=lapply(reord, function(ind){return(vert.ints.noreps[[ind]])})
}
get_genes_only <- function(featnames) {
	#print(head(featnames))
	#split.mat=get_split_mat(featnames)
	#return(split.mat[1,])
	list.elts=strsplit(as.character(featnames), "_")
	vec.of.genes=sapply(list.elts, function(x) {return(x[1])})
	names(vec.of.genes)=featnames
	return(vec.of.genes)
}

get_chrstates_only <- function(featnames) {
	#split.mat=get_split_mat(featnames)
	#return(split.mat[2,])
	list.elts=strsplit(as.character(featnames), "_")
	vec.of.states=sapply(list.elts, function(x) {return(x[2])})
	names(vec.of.states)=featnames
	return(vec.of.states)
}



get_split_mat <- function(featnames) {
	split.mat=sapply(featnames, function(x) {return(strsplit(x, "_")[[1]])})
	return(split.mat)
}

get.ranks.with.ties.for.nas <- function(vec) {
	index.rank=rank(vec, na.last=TRUE)
	na.inds=which(is.na(vec))
	last.rank=length(vec)
	index.rank[na.inds]=last.rank
	return(index.rank)
}
get.state.ecdf <- function(ranks, totalfeatures, rawnum=FALSE){
	xcoords=ranks
	if(rawnum) {
		ycoords=0:length(ranks)
	} else {
		ycoords=(0:length(ranks))/totalfeatures
	}
    return(stepfun(xcoords, ycoords))
}


plot_vector_as_hm<-function(outfile, data, rowlabels, collabels) { 
	pdf(outfile) 
	##image plots transpose so we transpose first
	par( mar = par( "mar" ) + c( 2, 4, 0, 0 ) ) 
	image(t(data), xaxt= "n", yaxt= "n" ) 
	axis( 1, at=seq(0,1,length.out=ncol( data ) ), labels=collabels, las= 2, cex.axis=.5) 
	axis( 2, at=seq(0,1,length.out=nrow( data ) ), labels=rowlabels, las= 2, cex.axis=.7 ) 
	dev.off() 
} 
