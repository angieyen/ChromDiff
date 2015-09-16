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
source("mkplots.R", chdir=T)
metric="perc"
testone=FALSE
covar.corr.off=FALSE

if(testone) {
	property="type"
	a.label="CellLine"
	b.label="PrimaryCulture"
	test="wilcox"
	test.correction="fdr"
	#window=1
	model="core_gencode_v10"
	covar.corr.off=FALSE
	metadatafile="data/final_celltype_metadata.txt"	
	#genefile=paste0("data/gencode_genes_", window, "kb_up.txt")
	genefile="data/gencode_genes_full.txt"
	covariate_mat_file="data/cov.mat.txt"
	map_covariates_file="data/map_vars_covariates.txt"
	expfile="data/57epigenomes.RPKM.pc"
	state_annotations_file="data/core_annotation.txt"	
	#generegions_label=paste0("gencode_v10_", window, "kb")
	generegions_label="gencode_v10"
	random=FALSE
	sampled=FALSE
	expronly=FALSE
	covar.corr.off=FALSE
}else{
	args<-commandArgs(TRUE)
	property=args[1]
	a.label=args[2]
	b.label=args[3]
	test=args[4]
	test.correction=args[5]
	model=args[6]
	metadatafile=args[7]
	genefile=args[8]
	covariate_mat_file=args[9]
	map_covariates_file=args[10]
	expfile=args[11]
	state_annotations_file=args[12]
	generegions_label=args[13]
}
## for randomization - use a.label="rand.${a.option}.NUM", b.label="rand.${b.option}.NUM", where NUM is unique ID
## for sampling - use a.label="sampled.${a.option}.NUMA", b.label="sampled.${b.option}.NUMB", where NUMA, NUMB, are the number of samples for group A and B
## for expronly - use a.label="expr.${a.option}", b.option="expr.${b.option}"

set_variables(model, metadatafile,  genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label)
a.option=parse_each_label(a.label)$option
b.option=parse_each_label(b.label)$option
pair.label=get_pair_label(a.label, b.label)
celltypes.list=get_celltypes(property, pair.label)
a.vec=celltypes.list$a.vec
b.vec=celltypes.list$b.vec

comp.groups(property, a.vec, b.vec, a.label, b.label, metric, test, test.correction, model, no_covariate_correction=covar.corr.off)

make.plots(property, a.vec, b.vec, a.label, b.label, metric, test, test.correction, model, no_covariate_correction=covar.corr.off)
