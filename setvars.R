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

assign("HOMEDIR", paste0(getwd(), "/"), envir=.GlobalEnv)



#### SET VARIABLES AND PARAMETERS #######
set_variables=function(curr_label, metadatafile, genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label) {
	assign("metadatafile", metadatafile, envir=.GlobalEnv)
	assign("curr_label", curr_label, envir=.GlobalEnv)
	assign("genefile", genefile, envir=.GlobalEnv)
	assign("covariate_mat_file", covariate_mat_file, envir=.GlobalEnv)
	assign("map_covariates_file", map_covariates_file, envir=.GlobalEnv)
	assign("expfile", expfile, envir=.GlobalEnv)
	assign("state_annotations_file", state_annotations_file, envir=.GlobalEnv)
	assign("generegions_label", generegions_label, envir=.GlobalEnv)
}

ensembl_dataset="hsapiens_gene_ensembl"

source("funcs2.R", chdir=T)

get_metadata=function() {
	rdatadir=paste("rdata/", curr_label, "/", sep="")
	metdata_rdata=paste(rdatadir,  "/celltype_metadata.rdata", sep="")
	if(file.exists(metdata_rdata) && file.info(metdata_rdata)$size>0) {
		load(metdata_rdata)
	} else {
		metdata=read.table(metadatafile, sep="\t", na.strings=c("NA", ""), header=TRUE, row.names=1, comment.char="")
		metdata=apply(metdata, 2, toupper)
		colnames(metdata)=toupper(colnames(metdata))
		kept_metdata=apply(metdata, 2, toupper)
		dir.create(rdatadir, recursive=TRUE, showWarnings=FALSE)
		save(kept_metdata, file=metdata_rdata)
	}
	return(kept_metdata)
}
get_validcelltypesfile=function(model, metric, property) {
	validcelltypes.rdatafile=paste0("rdata/", model, "/", metric,  "/validcelltypes.rdata")
	return(validcelltypes.rdatafile)
}
get_all_genes <- function() {
	genetable=get_geneinfo()
	allgenes=as.vector(genetable[, "id"])
	return(allgenes)
}

parse_labels <- function(pair.string) {
	groupelts=strsplit(pair.string, fixed=TRUE, split=".")[[1]]
	if((length(groupelts)%%2)!=0) {
		stop("Error: we expect the concatenation of the two labels to be of even length")
	} else {
		half=length(groupelts)/2
		label1=paste(groupelts[1:half], collapse=".")
		label2=paste(groupelts[(half+1):length(groupelts)], collapse=".")
	}
	return(list(label1=label1, label2=label2))
}


parse_each_label = function(label.string) {
	elts=strsplit(label.string, split=".", fixed=TRUE)[[1]]
	expr=FALSE
	random=FALSE
	sampled=FALSE
	num=NULL
	if((length(elts)==1)) {
		option=label.string
	} else {
		prefix=elts[1]
		if(!(prefix=="expr" || prefix=="rand" || prefix=="sampled")) {
			option=label.string
		}
		if(length(elts)==2 && prefix=="expr") {
			expr=TRUE
			option=elts[2]
		} else if (length(elts)==3) {
			if(prefix=="rand") {
				random=TRUE
			} else if(prefix=="sampled") {
				sampled=TRUE
			} else {
				stop(paste("Cannot parse ", pair.string))
			}
			option=elts[2]

			num=as.numeric(elts[3])
		}
	}
	return(list(option=option, expr=expr, random=random, sampled=sampled, num=num))
}

get_property <- function(pair.string) {
	type1str=parse_labels(pair.string)$label1
	type2str=parse_labels(pair.string)$label2
	
	list1=parse_each_label(type1str)
	list2=parse_each_label(type2str)
	if(list1$expr!=list2$expr || list1$random!=list2$random || list1$sampled!=list2$sampled) {
		stop(paste("Unable to parse", pair.string))	
	}

	type1str=list1$option
	type2str=list2$option
	kept_metdata=get_metadata()
	possible.types=colnames(kept_metdata)
	type1=tolower(names(which(apply(kept_metdata, 2, function(col) { return(toupper(type1str) %in% col)}))))
	type2=tolower(names(which(apply(kept_metdata, 2, function(col) { return(toupper(type2str) %in% col)}))))
	if(length(type1)==1 && (type1 %in% type2)) {
		type2=c(type1)
	} else if(length(type2)==1 && (type2 %in% type1)) {
		type1=c(type2)
	} else if(length(type1)>1 && length(type2)>1) {
		## remove specialgi 
		type1=type1[which(!(type1 %in% c("specialgi")))]
		type2=type2[which(!(type2 %in% c("specialgi")))]
	}

	if(length(type1)==1 && length(type2)==1 && type1==type2) {
		return(type1)
	} else {
		#save(list=ls(all=TRUE), file="tmp.rdata")
		stop(paste("No type found for", pair.string))
	}
}

get_ep_names = function(epids) {
	kept_metdata=get_metadata()
	ep_names=as.vector(kept_metdata[epids, "NAME"])
	return(ep_names)
}

get_ep_colors = function(epids) {
	kept_metdata=get_metadata()
	hex.colors=as.vector(kept_metdata[epids, "COLOR"])
	hex.colors[which(is.na(hex.colors))]="#000000"
	color_names=sapply(hex.colors, hex2col)
	return(color_names)
}

perform.wilcoxtest <- function(matrix.a, matrix.b) { 
    if(dim(matrix.a)[1]!=dim(matrix.b)[1]) { 
        stop("ERROR: data matrices must have same number of features (rows)") 
    } 
    numfeats=dim(matrix.a)[1] 
    p.vals=rep(NA, times=numfeats) 
    names(p.vals)=rownames(matrix.a) 
    for(feat.ind in 1:numfeats) { 
        if (feat.ind%%10000==0) { 
            print(paste("Processing feature ", feat.ind, "...", sep="")) 
        } 
        v1=matrix.a[feat.ind,] 
        v2=matrix.b[feat.ind,] 
        result=try(wilcox.test(matrix.a[feat.ind,], matrix.b[feat.ind,], alternative="two.sided", paired=FALSE, mu=0, na.action=na.omit), silent=TRUE) 
        if(!is(result,"try-error") && !is.na(result$p.val)) { 
            ## fill in p-values for successful tests 
            p.vals[feat.ind]=result$p.val 
        } 
    } 
    return(p.vals) 
} 
perform.ttest <- function(matrix.a, matrix.b) {
    if(dim(matrix.a)[1]!=dim(matrix.b)[1]) {
        stop("ERROR: data matrices must have same number of features (rows)")
    }
    numfeats=dim(matrix.a)[1]
    p.vals=rep(NA, times=numfeats)
    names(p.vals)=rownames(matrix.a)
    for(feat.ind in 1:numfeats) {
        if (feat.ind%%10000==0) {
            print(paste("Processing feature ", feat.ind, "...", sep=""))
    #       print(sort( sapply(ls(),function(x){object.size(get(x))})));
    #       print(gc())
    #       print(mem())
        }
        result=try(t.test(matrix.a[feat.ind,], matrix.b[feat.ind,], alternative="two.sided", paired=FALSE, mu=0, na.action=na.omit), silent=TRUE)
        if(!is(result,"try-error") && !is.na(result$p.val)) {
            ## append test statistic
            p.vals[feat.ind]= result$p.val
        }
    }
    return(p.vals)
}

perform.ftest <- function(matrix.a, matrix.b) {
    if(dim(matrix.a)[1]!=dim(matrix.b)[1]) {
        stop("ERROR: data matrices must have same number of features (rows)")
    }
    numfeats=dim(matrix.a)[1]
    p.vals=rep(NA, times=numfeats)
    names(p.vals)=rownames(matrix.a)
    for(feat.ind in 1:numfeats) {
        if (feat.ind%%10000==0) {
            print(paste("Processing feature ", feat.ind, "...", sep=""))
        }
        result=try(var.test(matrix.a[feat.ind,], matrix.b[feat.ind,], alternative="two.sided", na.action=na.omit), silent=TRUE)
        if(!is(result,"try-error") && !is.na(result$p.val)) {
            ## append test statistic
            p.vals[feat.ind]=result$p.value
        }
    }
    return(p.vals)
}

perform.bonferroni <- function(p.vals) {
    p.vals[which(is.na(p.vals))] = 1
    return(sort(p.vals*length(p.vals), decreasing=FALSE))
}
perform.fdr <- function(p.vals) {
    p.vals[which(is.na(p.vals))] = 1
    corrected.pvals=sort(p.adjust(p.vals, method="BH"))
    return(corrected.pvals)
}
perform.BY <- function(p.vals) {
    p.vals[which(is.na(p.vals))] = 1
    corrected.pvals=sort(p.adjust(p.vals, method="BY"))
    return(corrected.pvals)
}

get_featnames=function(model, metric, homedir="") {
    feat.dir=get_featdir(model)
    curr.feat.file=paste(homedir, feat.dir, get_featfile(metric), sep="")
    featnames =scan(curr.feat.file, what=character())
    return(featnames)
}


get_featdir<- function(model) { return(paste("featurenames/", model, "/", sep=""))}
get_rdatadir=function(model) {
	rdatadir=paste("rdata/", model, "/", sep="")
	return(rdatadir)
}

## variables for metrics
## elements correspond to percentage, delta, percentage windows, delta percentage windows
metric.labels=c("perc", "deltas", "windows", "deltawindows")
suffixes=c("_numsonly", "_deltas", "_window_perc_numsonly", "_windowdeltas")
feat.files=c("featurenames.txt", "featurenames.txt", "window_featnames.txt", "window_featnames.txt")
metric.subdirs=c("perc/", "deltas/", "windows/", "deltawindows/")
names(suffixes)=metric.labels
names(metric.subdirs)=metric.labels
names(feat.files)=metric.labels

## variables for test
test.labels=c("ttest", "wilcox", "ftest")
test.subdirs=c("ttest/", "wilcox/", "ftest/")
test.funcs=c(perform.ttest, perform.wilcoxtest, perform.ftest)
names(test.subdirs)=test.labels
names(test.funcs)=test.labels

## variables for hypothesis correction
corr.labels=c("bonferroni", "fdr", "BY")
corr.subdirs=c("bonferroni/", "fdr/", "BY/")
names(corr.subdirs)=corr.labels
corr.funcs=c(perform.bonferroni, perform.fdr, perform.BY)
names(corr.funcs)=corr.labels

check_input=function(metrics, tests, test.corrections) {
    ## match input metrics to metric.labels
    for(metric in metrics) {
        if(!(metric %in% metric.labels)) {
            stop(paste("metric ", metric, " is not valid option of ", toString(metric.labels), sep=""))
        }
    }
    ## match input test to test.labels
    for(test in tests) {
        if(!(test %in% test.labels)) {
            stop(paste("test ", test, " is not valid option of ", toString(test.labels), sep=""))
        }
    }
    ## match input correction to test.corrections
    for(correction in test.corrections) {
        if(!(correction %in% corr.labels)) {
            stop(paste("correction ", correction, " is not valid option of ", toString(corr.labels), sep=""))
        }
    }
}
get_corr_func=function(correction) {
	return(corr.funcs[[correction]])
}
get_test_subdir=function(test) {
	return(test.subdirs[test])
}
get_test_func=function(test) {
	return(test.funcs[[test]])
}

get_suffix <- function(metric) {
	return(suffixes[metric])
}
get_prefix <- function(metric, model) {
	prefix1=paste("perc/", model, "/numsonly/", sep="")
	prefix2=paste("deltas/", model, "/", sep="")
	prefix3=paste("perc/", model, "/", sep="")
	prefix4=paste("deltas/", model, "/", sep="")
	prefixes=c(prefix1, prefix2, prefix3, prefix4)
	names(prefixes)=metric.labels
	return(prefixes[metric])
}
get_metric_subdir <- function(metric) {
	return(metric.subdirs[metric])
}

get_corr_subdir=function(correction) {
	return(corr.subdirs[correction])
}


get_featfile <- function(metric) {
	return(feat.files[metric])
}

get_covariates_to_correct = function() {
	mapping=read.table(map_covariates_file, row.names=1, sep="\t")
	covariates_to_correct=rownames(mapping)
	return(covariates_to_correct)
}

get_all_plottypes = function() {
	return(c("usedend", "assocstate", "domstate", "combinations"))
}

vec.match =function(vec1, vec2) {
    new.v1=sort(vec1)
	new.v2=sort(vec2)
	return(new.v1==new.v2)
}

get_plottype.labels = function() {
    return(c("usedend", "assocstate", "domstate", "combinations"))
}

get_plottype.suffixes = function() {
    arr=c()
    arr["assocstate"]="_assocstate"
    arr["usedend"]="_dend"
    arr["domstate"]="_domstate"
    arr["combinations"]="_combs"
    if(length(arr)!=length(get_all_plottypes()) || !vec.match(names(arr),get_all_plottypes())) {
        stop("Check get_plottype.labels() and get_plottype.suffixes(): Plottype labels and plottype suffixes do not match")
    }
    return(arr)
}

get_plottype_suffix=function(plottype) {
    if (!(plottype %in% get_all_plottypes())) {
        stop(paste("plottype ", plottype, " not recognized.", sep=""))
    }
    return(get_plottype.suffixes()[plottype])
}

get_plottype_bools=function(plottype) {
	reorder_by_assocstate=FALSE
	useDend=FALSE
	reorder_by_domstate=FALSE
	reorder_by_combinations=FALSE
	if(plottype=="assocstate") {
		reorder_by_assocstate=TRUE
	} else if (plottype=="usedend"){
		useDend=TRUE
	} else if (plottype=="domstate") {
		reorder_by_domstate=TRUE
	} else if(plottype=="combinations") {
		reorder_by_combinations=TRUE
	} else {
		stop(paste("plottype ", plottype, " not recognized.", sep=""))
	}
	return(list(reorder_by_assocstate, useDend, reorder_by_domstate, reorder_by_combinations))
}

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
	args <- as.list(match.call())
	args <- args[-c(1:2,length(args))]
	length(value) <- length(args)
	for(i in seq(along=args)) {
		a <- args[[i]]
		if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
	}
	x
}


get_metadata_properties=function() { 
	kept_metdata=get_metadata()
	return(colnames(kept_metdata))
}

##### get matrix of covariates ####
get_covariate_mat <- function(dependentvariable) {
	fixed.cov.mat=get_full_covariate_mat()

	dependentvariable=toupper(dependentvariable)
	covariate_columns_to_exclude=get_columns_for_variable(dependentvariable)
	new_colnames=colnames(fixed.cov.mat)[which(!(colnames(fixed.cov.mat) %in% covariate_columns_to_exclude ))]
	final.cov.mat=fixed.cov.mat[,new_colnames]
	return(final.cov.mat)
}

get_columns_for_variable<- function(dependentvariable) {
	if(file.exists(map_covariates_file)) {
		map.matrix=as.matrix(read.table(map_covariates_file, sep="\t"))
	} else {
		stop(paste("Mapping of metadata matrix to covariate matrix does not exist as specified by variable map_covariates_file: ", map_covariates_file))
	}
	if(dependentvariable %in% rownames(map.matrix)) {
		assoc.covar.cols=colnames(map.matrix)[which(map.matrix[dependentvariable,])]
	} else {
		assoc.covar.cols=c()
	}
	return(assoc.covar.cols)
}

get_full_covariate_mat <- function() {
	if(file.exists(covariate_mat_file)) {
		fixed.cov.mat=read.table(covariate_mat_file, sep="\t")
	} else {
		stop(paste("Covariate matrix does not exist as specified by variable covariate_mat_file: ", covariate_mat_file))
	}
	return(fixed.cov.mat)
}

get_geneinfo=function(homedir="") {
	geneinfo=read.table(paste0(homedir, genefile))
    colnames(geneinfo)=c("chr", "start", "stop", "strand", "id", "symbol")
    rownames(geneinfo)=geneinfo[, "id"]
	return(geneinfo)
}

get_genenames <- function(geneids, homedir="") {
	geneinfo=get_geneinfo(homedir)
	results=as.vector(geneinfo[geneids, "symbol"])
	return(results)
}
get_genesizes <- function(geneids, homedir="") {
	geneinfo=get_geneinfo(homedir)
	starts=geneinfo[geneids,"start" ]
	ends=geneinfo[geneids,"stop" ]
	sizes=ends-starts
	return(sizes)
}
get_ensembl.entrez.map = function(ensemblids) {
	library("biomaRt")
	## choose database
	ensemblMart=useMart("ensembl", dataset=ensembl_dataset)
	## choose inputtype from listFilters(ensemblMart) and outputtype from listAttributes(ensemblMart)
	inputtype="ensembl_gene_id"
	inputvalues=ensemblids
	outputtype=c("ensembl_gene_id", "entrezgene")
	result=getBM(attributes=outputtype, filters=inputtype, values=inputvalues, mart=ensemblMart)
	return(result)
}

get_entrezids <- function(ensemblids) {
	return(get_ensembl.entrez.map(ensemblids)[,"entrezgene"])
}

get_ensembl.ids = function(gencodeids) {
    ensembl.genes=sapply(gencodeids, function(x) {elts=strsplit(x, split=".", fixed=TRUE); return(elts[[1]][1])})
	return(ensembl.genes)
}

#special_optlabels=list("CELLLINEALL"=toupper(c("CellLineDerived", "CellLine", "CellLine_Cancer")), "CELLLINEANDDERIVED"=toupper(c("CellLineDerived", "CellLine")))

get_celltypes_helper <- function(property, options, expronly=FALSE) {
	kept_metdata=get_metadata()
	property=toupper(property)
	options=toupper(options)
	## make sure property is valid
	if (!(property %in% get_metadata_properties())) {
		stop(paste("property ", property, " is not in list of available options: ", toString(get_metadata_properties())))
	} else {
		## make sure all options are valid
		valid_options=unique(kept_metdata[,property])
		optionstoadd=c()
		standardopts=c()
		for (option in options) {
			if (!(option %in% valid_options)) {
				stop(paste("option ", option, " has no matches for property ", property, sep=""))
			} else {
				standardopts=append(standardopts, option)
			}
		}
	}

	alloptions=unique(append(standardopts, optionstoadd))
	celltypes=rownames(kept_metdata)[which(kept_metdata[, property] %in% alloptions)]

	epids_with_expr=rownames(get_all_logexp())
	if(expronly) { #filter only to epigenomes with expression data
		celltypes=celltypes[which(celltypes %in% epids_with_expr)]
	}

	return(celltypes) 
}

get_valid_celltypes=function(metric, property, labels.str) {
	valid_celltypesfile=get_validcelltypesfile(curr_label, metric, property)
	if(file.exists(valid_celltypesfile) && file.info(valid_celltypesfile)$size>0) {
		load(valid_celltypesfile)
	} else {
		stop("Can not access valid celltypes...")
	}
	result=get_celltypes(property, labels.str)
	kept.list1=result$a.vec[which(result$a.vec %in% validcelltypes)]
	kept.list2=result$b.vec[which(result$b.vec %in% validcelltypes)]
	return(list(a.vec=kept.list1, b.vec=kept.list2))
}


 
get_celltypes=function(property, labels.str) {
    labels.list=parse_labels(labels.str)
    a.label=labels.list$label1
    b.label=labels.list$label2
    parsed.label1=parse_each_label(a.label)
    parsed.label2=parse_each_label(b.label)
    a.option=parsed.label1$option
    b.option=parsed.label2$option
	if(parsed.label1$expr!=parsed.label2$expr || parsed.label1$random!=parsed.label2$random || parsed.label1$sampled!=parsed.label2$sampled) {
    	stop(paste("Unable to parse", pair.string))
    }
	expronly=parsed.label1$expr
	random=parsed.label1$random
	sampled=parsed.label1$sampled
	num_a=parsed.label1$num
	num_b=parsed.label2$num
	
	if(sampled && random) {
    	stop("Can not have both sampling and randomization turned on.")
	}

	a.vec=get_celltypes_helper(property, c(a.option), expronly)
	b.vec=get_celltypes_helper(property, c(b.option), expronly)
	if(sampled) {
    	sampled.dir=paste("sampled_assignments/", curr_label, "/", sep="")
    	sampled.file=paste(sampled.dir, a.label, ".", b.label, ".rdata", sep="")
    	if(file.exists(sampled.file)) {
			load(sampled.file)
    	} else {
        	a.elts=strsplit(a.label, split=".", fixed=TRUE)[[1]]
        	b.elts=strsplit(b.label, ".", fixed=TRUE)[[1]]
        	check=function(vec) {
            	if(length(vec)!=3 || vec[1]!="sampled" || !is.integer(as.integer(vec[3]))) {
                	stop(paste("Unable to sample based on label: ", vec,". Please use format sampled.OPTION.NUM to sample group OPTION down to NUM epigenomes",  sep=""))
            	}
        	}
        	check(a.elts)
        	check(b.elts)

        	num_a=as.integer(a.elts[3])
        	num_b=as.integer(b.elts[3])
     	 
        	sampled.a.vec=sample(a.vec, size=num_a, replace=FALSE)
        	sampled.b.vec=sample(b.vec, size=num_b, replace=FALSE)
        	dir.create(sampled.dir, recursive=TRUE, showWarnings=FALSE)	
        	save(list=c("sampled.a.vec", "sampled.b.vec"), file=sampled.file)
    	}
    	a.vec=sampled.a.vec
    	b.vec=sampled.b.vec
	} else if(random) {
    	rand.dir=paste("rand_assignments/", curr_label, "/", sep="")
    	rand.assign.file=paste(rand.dir, a.label, ".", b.label, ".rdata", sep="")
    	if(file.exists(rand.assign.file)) {
        	load(rand.assign.file)
    	} else {
        	all=c(a.vec, b.vec)
        	rand.a=sample(x=all, size=length(a.vec), replace=FALSE)
        	b.inds=which(!(all %in% rand.a))
        	rand.b=all[b.inds]
        	#print(rand.a)
        	#print(rand.b)
        	a.vec=rand.a
        	b.vec=rand.b
        	dir.create(rand.dir, recursive=TRUE, showWarnings=FALSE)	
        	save(list=c("a.vec", "b.vec"), file=rand.assign.file)
    	}
	}
	return(list(a.vec=a.vec, b.vec=b.vec))
}

get_all_logexp=function() {
	expdir=paste("expression/", curr_label, "/", sep="")
    rdatafile=paste(expdir,"expressionmat.rdata", sep="")
    if(!file.exists(rdatafile)) {
        #print("Reading exp data...")
        ## in orig file: columns are celltypes, ensembl gene ids are rows
        origmat=read.table(expfile, row.names=1, header=TRUE)
        expmat=t(origmat)
        dir.create(expdir, recursive=TRUE, showWarnings=FALSE)
        save(expmat, file=rdatafile)
    } else {            
        #print("Loading exp data...")
        load(rdatafile) 
    }                       
    logexp.mat=log(expmat+1)
    return(logexp.mat)
}

get_exp_resids=function(property, residsdir) {
	## check if dependent variable is one of the covariates we correct for
	if(toupper(property) %in% get_covariates_to_correct()) {
    	varsuffix=paste(".skip", property, sep="")
	} else {
    	varsuffix=".skipnone"
	}
	residsfile=paste0(residsdir, "residmat", varsuffix, ".rdata")
	cov.mat = get_covariate_mat(property)

	## get corrected and uncorrected data, as well as random data and data for nonsig genes
	logexp.mat=get_all_logexp()
	geneswithinfo=colnames(logexp.mat)[which(apply(logexp.mat, 2, function(x) {return(!(all(is.na(x))))}))]
	epswithexpdata=rownames(logexp.mat)[which(apply(logexp.mat, 1, function(x) {return(!(all(is.na(x))))}))]
	## throw out eps without covariate data (they were blacklisted data out for other reasons)
	kept_epswithexpdata=epswithexpdata[which(epswithexpdata %in% rownames(cov.mat))]
	logexp.mat=logexp.mat[kept_epswithexpdata, geneswithinfo]
	## correct matrix based on covariate mat
	if(file.exists(residsfile) && file.info(residsfile)$size>0) {
    	load(residsfile)
	} else {
    	## filter cov matrix for eps with exp data
    	fixed.cov.mat=cov.mat[kept_epswithexpdata, ]
    	resids.mat=get_residuals_from_mats(fixed.cov.mat, logexp.mat, method="lin")
    	dir.create(residsdir, recursive=TRUE, showWarnings=FALSE)
    	save(resids.mat, file=residsfile)
	}
	return(resids.mat)
}

get_random_celltypes <- function(ngroupa, ngroupb) {
	celltypes=get_all_celltypes()
	inds=sample(1:length(celltypes), ngroupa+ngroupb, replace=F)
	celltypesa=celltypes[inds[1:ngroupa]]
	celltypesb=celltypes[inds[(ngroupa+1):(ngroupa+ngroupb)]]
	return(list(celltypesa, celltypesb))
}
get_all_celltypes <- function() {
	kept_metdata=get_metadata()
	celltypes=rownames(kept_metdata)
	return(celltypes)
}

#get_valid_celltypes <- function(model) {
#	valid.celltypes=c()
#	celltypes.list=get_all_celltypes()
#	for (ind in 1:length(celltypes.list)) {
#		datafile=paste("deltas/", model, "/",  celltypes.list[ind], "_deltas.txt", sep="")
#		if(file.exists(datafile)) {
#			valid.celltypes=append(valid.celltypes, celltypes.list[ind])
#		}
#	}
#	return(valid.celltypes)
#}

get_gene_symbols=function(geneids) {
    allgenes=read.table(genefile)
    result=as.vector(allgenes[geneids, "symbol"])
    names(result)=geneids
    return(result)
}

mapStateNumToName <- function(statevec) {
	annot.info=get_annot.info()
	statenames=annot.info[, "MNEMONIC"]
	names(statenames)=rownames(annot.info)
	return(statenames[statevec])
}
mapStateNameToNum <- function(statenamevec) {
	annot.info=get_annot.info()
	statenums=rownames(annot.info)
	statenames=annot.info[, "MNEMONIC"]
	names(statenums)=statenames
	return(statenums[statenamevec])
}

get_all_chrstates = function() {
	annot.info=get_annot.info()
	chrstates=as.numeric(rownames(annot.info))
	return(chrstates)
}
get_chrstate_colors <- function(chrstate.vec) {
	default="mistyrose1"
	allcolors=get_state_colors()
	chromatin.labels=rep(default, length(chrstate.vec)) 
	names(chromatin.labels)=chrstate.vec
	for(currstate in unique(as.character(chrstate.vec))) {
		if(!is.na(currstate)) {
			inds=which(names(chromatin.labels)==currstate)
			chromatin.labels[inds]=allcolors[currstate]
		}
	}
	if(length(which(chromatin.labels==default))>0) {
		warning("chromatin state vector did not map correctly to chromatin state colors")
	}
	return(chromatin.labels)
}

get_annot.info <- function() {
    annot.info=read.table(state_annotations_file, header=TRUE, row.names=1, sep="\t")
    return(annot.info)
}

get_state_colors <- function(repWhite=FALSE){
    annot.info=get_annot.info()
    codes=as.character(annot.info[, "COLOR.CODE"])
    names(codes)=rownames(annot.info)
    color.code.matrix=sapply(as.character(codes), function(x) {unlist(strsplit(x, ","))})
    colors=rgb(color.code.matrix[1,], color.code.matrix[2,], color.code.matrix[3,], maxColorValue=255)
    names(colors)=rownames(annot.info)
    ## replace white with black
    if(repWhite) {
        colors[which(colors=="#FFFFFF")]="#000000"
    }
    return(colors)
}

get_plotdir=function(model, a.label, b.label) {
    pair.label=get_pair_label(a.label, b.label)
    base.plotdir=paste("plots/", model,"/", sep="")
    plotdir=paste(base.plotdir, pair.label, "/", sep="")
    return(plotdir)
}

get_full_rdatadir=function(model, test, metric) {
	rdatadir=get_rdatadir(model)
	test.subdir=get_test_subdir(test)
	metric.subdir=get_metric_subdir(metric)
	full.rdatadir=paste(rdatadir, metric.subdir, test.subdir, sep="")
	return(full.rdatadir)
}
get_plot_subdir=function(test, metric, correction) {
	test.subdir=get_test_subdir(test)
	metric.subdir=get_metric_subdir(metric)
	corr.subdir=get_corr_subdir(correction)
   	subdir=paste(metric.subdir, test.subdir, corr.subdir, sep="")
	return(subdir)
}

get_full_plotdir=function(model, a.label, b.label, test, metric, correction) {
	plotdir=get_plotdir(model, a.label, b.label)
	plot_subdir=get_plot_subdir(test, metric, correction)
	full_plotdir=paste(plotdir, plot_subdir, sep="")
	return(full_plotdir)
}

get_pair_label=function(a.label, b.label) {
	return(paste(a.label, b.label, sep="."))
}

get_identifying_prefix=function(metric, test, correction, property, a.label, b.label) {
	pair.label=get_pair_label(a.label, b.label)
	prefix=paste(metric, test, correction, property, pair.label, sep=".")
	return(prefix) 
}

get_geneorder_file = function(plotdir, suffix, plottypesuffix) { 
    geneorderfile=paste(plotdir, "/sig_maj_geneorder", suffix, plottypesuffix, ".txt", sep="") 
	return(geneorderfile) 
}

get_expdata_file=function(plotdir, plottypesuffix) {
	expdatafile=paste(plotdir, "exp.mat.to.plot", plottypesuffix, ".Rdata", sep="")
	return(expdatafile)
}

get_pval_file=function(model, metric, test, correction, property, a.label, b.label) {
	prefix=get_identifying_prefix(metric, test, correction, property, a.label, b.label)
	allpval.file=paste("all_pvals/", model, "/",prefix, ".txt", sep="")
	return(allpval.file)
}
get_sigpval_file=function(model, metric, test, correction, property, a.label, b.label) {
	prefix=get_identifying_prefix(metric, test, correction, property, a.label, b.label)
	sigpval.file=paste("sig_pvals/", model, "/", prefix, ".txt", sep="")
	return(sigpval.file)
}
