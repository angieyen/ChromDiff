source("setvars.R", chdir=T)
source("funcs.R", chdir=T)

## function for drawing clusters
final.drawing = function(horiz.ints, vert.ints, linewidth, minclustsize, total_rows) {
    #print(horiz.ints)
    #print(vert.ints)

    # draw horizontal line
    abline(h=horiz.ints,  lwd=linewidth)
    ## draw rectangles for groups that meet the cutoff
    drawRects(vert.ints, linewidth, minclustsize, total_rows)
}

drawRects <- function(vert.ints, linewidth, minclustsize, nrows) {
    ybottom=0.5
    ytop=nrows+.5
    if(length(vert.ints)>0) {
        for(ind in 1:length(vert.ints)) {
            currclust=vert.ints[[ind]]
            if((currclust[2]-currclust[1])>=minclustsize) {
                rect(currclust[1]-.5, ybottom, currclust[2]+.5, ytop, lwd=linewidth)
            }
        }
    }
}

get_vert.ints_dend = function(plotdir, mat.to.plot, noreps=FALSE) {
	dendfile=paste(plotdir,"sigfeats_feats_dend_matched.rdata", sep="")
	d1file=paste(plotdir, "sigfeats_feats_d1_matched.rdata", sep="")
	d3file=paste(plotdir, "sigfeats_feats_d3_matched.rdata", sep="")
    groupings.file=paste(plotdir,"sigfeats_feats_colname_groupings_matched.rdata", sep="")
    ## if no repeat genes
    if(!file.exists(d1file)) {
        if(!file.exists(dendfile)) {error(paste("Neither file ", d1file, " nor ", dendfile, " exist.", sep=""))}
        load(dendfile)
        vert.ints=get_vert_ints_height(colnames(mat.to.plot),dend, heightcutoff)
        vert.ints.noreps=vert.ints
        #d1offset=0
        #d2offset=0
    } else {
    ## if there were repeat genes
        load(groupings.file)
        sigfeats=colnames(mat.to.plot)
        totalsize=length(sigfeats)
        #print(sigfeats)
        #print(groupings)
        ## Use the saved dendrogram and pvclust info if possible
        if(file.exists(d1file)) {
            load(d1file)
            #print(d1file)
            d1offset=length(final.group1.featnames)
            if(length(final.group1.featnames)>1) {
                d1.vert.ints=get_vert_ints_height(colnames(mat.to.plot), d1, heightcutoff)
                if(d1offset!=length(labels(d1))) { warning("There should be the same length of elements in d1 dendrogram as final.group1.featnames")}
            } else {
                #print(length(final.group1.featnames))
                ## there should only be 1 feature in group 1 in this case
                if(d1offset!=1) {warning("There should only be 1 feature in group 1...")}
                d1.vert.ints=list(c(1,d1offset))
            }
        } else {
        ## use grouping information directly    
            warning("No saved dendrogram but using grouping information for group 1...")
            d1offset=length(which(groupings==1))
            if(d1offset>0) {
                d1.vert.ints=list(c(1,d1offset))
            }
        }

        ### get d2 ints based on gene names
        group2feats=names(groupings)[which(groupings==2)]
        d2offset=length(group2feats)
        if(d2offset>0) {
            group2genes=get_genes_only(group2feats)
            d2.vert.ints=list()
            ## figure out which genes are right before divisions (last genes in cluster)
            for(ind in 1:length(d1.vert.ints)) {
                dividing.genes = get_genes_only(sigfeats[d1.vert.ints[[ind]]])
                dividing.inds=match_last(dividing.genes, group2genes)
                d2.vert.ints[[ind]]=dividing.inds+d1offset
            }
        } else {
            warning("Group 2 should not have a 0 or negative offset")
        }

        d3offset=length(which(groupings==3))
        if(d3offset==0) {
            d3.vert.ints=list()
        } else {
            if(file.exists(d3file)) {
                load(d3file)
                if(d3offset!=length(final.group3.featnames)){warning("Group 3 offset is not consistent...") }
                if(length(final.group3.featnames)>1) {
                    d3.vert.ints=get_vert_ints_height(colnames(mat.to.plot), d3, heightcutoff)
                } else if(d3offset==1) {
                    if (length(final.group3.featnames)!=1){ warning("There should be only 1 feature in group 3...")}
                    d3.vert.ints=list(c((d1offset+d2offset+1),(d1offset+d2offset+1)))
                }
            } else {
                ## use grouping information directly    
                warning("No saved dendrogram but using grouping information for group 1...")
                d3.vert.ints=list(c(d1offset+d2offset+1,d1offset+d2offset+d3offset))
            }
        }
        if(d1offset + d2offset+d3offset != totalsize) {
            warning("Error with grouping sizes")
        }
        vert.ints=c(d1.vert.ints, d2.vert.ints, d3.vert.ints)
        d3.vert.ints.noreps=lapply(d3.vert.ints, function(x){return(x-d2offset)})
        vert.ints.noreps=c(d1.vert.ints, d3.vert.ints.noreps)
    }
    if(noreps) {
    	returnlist=list(vert.ints.noreps=vert.ints.noreps, mat.to.plot=mat.to.plot)
    } else {
    	returnlist= list(vert.ints=vert.ints, mat.to.plot=mat.to.plot)
    }
	return(returnlist)
}

get_vert.ints_assocstate = function(plotdir, mat.to.plot) {
    ## to order by associated chr state
    colstates = as.numeric(get_chrstates_only(colnames(mat.to.plot)))
    new_order=order(colstates)
    mat.to.plot=mat.to.plot[,new_order]
    first.inst=sapply(min(colstates):max(colstates), function(x) {which(colstates[new_order]==x)[1]})
    first.inst=first.inst[which(!is.na(first.inst))] 
    vert.ints=list()
    for(ind in 1:length(first.inst)) {
        if(ind<length(first.inst)) { 
            vert.ints[[ind]]=c(first.inst[ind], first.inst[ind+1]-1)
        } else {
            vert.ints[[ind]]=c(first.inst[ind], ncol(mat.to.plot))
        }   
    }   
    return(list(vert.ints=vert.ints, mat.to.plot=mat.to.plot))
}
get_vert.ints_domstate = function(plotdir, mat.to.plot) {
    ## rename matrix columns so that they are gene names instead of feature ids
    feat2genes=get_genes_only(colnames(mat.to.plot))
    names(feat2genes)=colnames(mat.to.plot)
    colnames(mat.to.plot)=feat2genes
    
    ## load in gene order from previous sig majority plotting
    suffix="_matched"
    plottypesuffix="_domstate"
    geneorderfile=paste(plotdir, "sig_maj_geneorder", suffix, plottypesuffix, ".txt", sep="")
    geneorder=as.vector(t(read.table(geneorderfile, header=TRUE)))
    ## reorder matrix based on previous sig majority plotting for domstate
    mat.to.plot=mat.to.plot[, geneorder]
    ## load previously calculated dendrogram
    majdendfile=paste(plotdir, "sig_maj_plot_dend", plottypesuffix, ".rdata", sep="")
    load(majdendfile) 
    ## get clusterings based on previous dendrogram
    vert.ints=reorder_get_vert_ints_height(colnames(mat.to.plot), dend, heightcutoff)
    
    ## convert the colnames back to features
    genes2feats=names(feat2genes)
    names(genes2feats)=feat2genes
    colnames(mat.to.plot)=genes2feats[colnames(mat.to.plot)]
	
	return(list(vert.ints=vert.ints, mat.to.plot=mat.to.plot))
}

get_vert.ints.file=function(plotdir, plottypesuffix) {
	vert.ints.file=paste(plotdir, "sigfeats_feats_clust_vert_ints", plottypesuffix, ".rdata", sep="") 
	return(vert.ints.file)
}
get_vert.ints_combinations = function(plotdir, mat.to.plot) {
    sigpvals.order.file=paste(plotdir, "sig.pvals.rdata", sep="")
    load(sigpvals.order.file)
    
    ## rename matrix columns so that they are gene names instead of feature ids
    feat2genes=get_genes_only(colnames(mat.to.plot))
    names(feat2genes)=colnames(mat.to.plot)
    colnames(mat.to.plot)=feat2genes
    
    ## reorder based on previous ordering
    mat.to.plot=mat.to.plot[, sigpvals.order]
	###FIX THIS###
	#Don't have clusterings because don't have dendrogram (only have ordering)
	#dend=sigpvals.dend
	#vert.ints=reorder_get_vert_ints_height(colnames(mat.to.plot), dend, heightcutoff)
	vert.ints=list()

    ## convert the colnames back to features
    genes2feats=names(feat2genes)
    names(genes2feats)=feat2genes
    colnames(mat.to.plot)=genes2feats[colnames(mat.to.plot)]

    return(list(vert.ints=vert.ints, mat.to.plot=mat.to.plot))
}

get_vert.ints = function(plottype.str, plotdir, mat.to.plot) {
	list[reorder_by_assocstate, reorder_by_dend, reorder_by_domstate, reorder_by_combinations] = get_plottype_bools(plottype.str)

	if(reorder_by_dend) {
		ls= get_vert.ints_dend(plotdir, mat.to.plot)
	} else if(reorder_by_assocstate){
		ls=get_vert.ints_assocstate(plotdir, mat.to.plot)
	} else if(reorder_by_domstate) {
		ls=get_vert.ints_domstate(plotdir, mat.to.plot)
	} else if(reorder_by_combinations) {
		ls= get_vert.ints_combinations(plotdir, mat.to.plot) 
	}
	return(ls)
}
