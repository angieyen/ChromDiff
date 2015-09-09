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

source("../setvars.R", chdir=T)
library("qvalue") 
gencodemappedfile=paste0(HOMEDIR, "genes/", generegions_label, "/gencode_mapped_symbols.bed")
## number of all human gene symbols taken from HUGO gene symbols on Feb 2014
get_msigdata = function(msigdbset, idtype){
    rdatafile=paste0(HOMEDIR,"msigdb/rdata/", msigdbset,".", idtype, ".rdata")
    if(file.exists(rdatafile)) {
        load(rdatafile)
    } else {
        msigfile=paste0(HOMEDIR, "msigdb/genesets/", msigdbset, ".all.v4.0.", idtype, ".gmt")
        numcols=count.fields(msigfile)
        msigdata=read.table(msigfile, fill=TRUE, col.names=1:max(numcols))
        save(msigdata, file=rdatafile)
    }
   ### FIX THIS HERE #### 
    ## remove any symbols not mapped from gencode
    if(idtype=="symbols") {
    	filt.msigdatafile=paste(HOMEDIR,"msigdb/rdata/", msigdbset,".", idtype, ".filtered.rdata", sep="")
    	if(file.exists(filt.msigdatafile)) {
    		load(filt.msigdatafile) 
    	} else {
			mappedsymbols=as.vector(t(read.table(gencodemappedfile)))
    		newmsigdata=apply(msigdata, 1, function(msigdbgroup) {
    			currgroupgenes=msigdbgroup[3:length(msigdbgroup)]
    			indstoconvert=which(!(currgroupgenes %in% mappedsymbols) && currgroupgenes!="")
    			currgroupgenes[indstoconvert]=""
    			return(c(msigdbgroup[1:2], currgroupgenes))
    		})
    		msigdata=t(newmsigdata)
    		save(msigdata, file=filt.msigdatafile)
    	}	
    }
    return(msigdata)
}

get_sig_genes=function(pval.table) {
	sig_feats=which(pval.table$pval<.05)
	sig_pval.table=pval.table[sig_feats,]
	sig_genes=get_genes_only(sig_pval.table$feat)
	return(sig_genes)
}

get_msigdb_resultsdir=function(homedir,msigdbset, idtype, model) {
    resultsdir=paste(homedir, "msigdb/results/", msigdbset, "/",idtype, "/", model, "/", sep="")
	return(resultsdir)
}
get_msigdb_resultsfile=function(homedir, msigdbset, idtype, model, metric, test, correction, property, a.label, b.label) {
    resultsdir=get_msigdb_resultsdir(homedir, msigdbset, idtype, model)
    prefix=get_identifying_prefix(metric, test, correction, property, a.label, b.label)
    outfile=paste(resultsdir, prefix,"_all.txt", sep="")
	return(outfile)
}	
get_msigdb_clustfile=function(homedir, msigdbset, idtype, model, metric, test, correction, property, a.label, b.label, clustcount) {
    resultsdir=get_msigdb_resultsdir(homedir, msigdbset, idtype, model)
    prefix=get_identifying_prefix(metric, test, correction, property, a.label, b.label)
	outfile=paste(resultsdir,prefix, "_subgroup", clustcount, ".txt", sep="")
	return(outfile)
}	

get_gene.universe=function(idtype) {
	if(idtype=="entrez") {
		## number taken from the Broad?
		# tmp: set # of genes to what broad said
		gene.universe=45956
	} else if (idtype=="symbols") {
		## only consider universe as symbols that are mapped from gencode
		gene.universe=as.vector(t(read.table(gencodemappedfile)))
	}
	return(gene.universe)
}

convert=function(gencode.genes, idtype, homedir="") {
	if(idtype=="entrez") {
		ensembl.genes=get_ensembl.ids(gencode.genes)
		final.genes=convert_ensembl2entrez(ensembl.genes)
	} else if (idtype=="symbols") {
		final.genes=get_genenames(gencode.genes, homedir)
	}
	else {
		stop(paste(idtype, "not recognized by convert function from gencode genes"))
	}
	return(final.genes)
}	

convert_ensembl2entrez=function(ensembl.genes) {
    ensembl.entrez.mat=get_ensembl.entrez.map(ensembl.genes)
    entrezcol="entrezgene"
    ensemblcol="ensembl_gene_id"
	final.genes=ensembl.entrez.mat[,entrezcol]
	names(final.genes)=ensembl.entrez.mat[,ensemblcol]
	return(final.genes)
}
remove.extra=function(genes) {
    tmp=unique(genes[which(!is.na(genes))])
    return(tmp[which(tmp!="")])
}

#get_enrichments.ensembl=function(msigdata, ensembl.genes) {
#    entrez.genes=convert_ensembl2entrez(ensembl.genes)
#    results=get_enrichments(msigdata, entrez.genes)
#    return(results)
#}

#get_enrichments.gencode=function(msigdata, gencode.genes) {
#    ensembl.genes=get_ensembl.ids(gencode.genes)
#   	return(get_enrichments.ensembl(msigdbset, ensembl.genes))
#}


get_enrichments <- function(msigdata, genelist, gene.universe) {
    genesetsizes=apply(msigdata, 1, function(msigdbgroup) {return(length(which(msigdbgroup[3:length(msigdbgroup)]!="")))})

    numoverlaps=apply(msigdata, 1, function(msigdbgroup) { msigdbgenes=msigdbgroup[3:length(msigdbgroup)]; numoverlap=length(which(genelist %in% msigdbgenes)); return(numoverlap) })
    genelistsize=length(genelist)
    ## calculate hypergeom p-vals
    pvals=lapply(1:nrow(msigdata), function(ind) { 
    	k=numoverlaps[ind]-1
    	genesetsize=genesetsizes[ind]
    	return(phyper(k, genesetsize, gene.universe-genesetsize, genelistsize, lower.tail=FALSE))})
    pvals=unlist(pvals)
    names(pvals)=msigdata[,1]
    ordering=order(pvals, decreasing=FALSE)
    sorted.pvals=pvals[ordering]
	print(paste("Number of unique gene ids used for enrichment:", genelistsize))
    ## Use Benjamini-Hochberg FDR correction
    #fdr.vals=p.adjust(sorted.pvals, method="BH")
   	## Use Storey's q-value
   	fdr.vals=qvalue(sorted.pvals)$qvalue
   	names(fdr.vals)=names(sorted.pvals)

	results=data.frame(geneset=names(fdr.vals), genesetsizes=as.numeric(genesetsizes[ordering]), numoverlaps=as.numeric(numoverlaps[ordering]), koverK=as.numeric(numoverlaps[ordering])/as.numeric(genesetsizes[ordering]), qval=as.numeric(fdr.vals), origpval=as.numeric(sorted.pvals))	
	#results=cbind(names(fdr.vals), as.numeric(genesetsizes[ordering]), as.numeric(numoverlaps[ordering]), as.numeric(numoverlaps[ordering])/as.numeric(genesetsizes[ordering]), as.numeric(fdr.vals), as.numeric(sorted.pvals))	
	row.names(results)=names(fdr.vals)
    colnames(results)=c("geneset", "genesetsize (K)", "numoverlaps (k)", "k/K", "qval", "orig_pval")
    keepinds=which(results[,"qval"]<.05)
    results=results[keepinds,]
    if(length(keepinds)>1) {
        tmporder=order(results[, 4], decreasing=TRUE)
        tmp=results[tmporder,]
        finalorder=order(tmp[, "qval"], decreasing=FALSE)
        results=tmp[finalorder,]
    }
    return(results)
}

