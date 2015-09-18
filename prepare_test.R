source("funcs.R")
model="core_gencode_v10"
metadatafile="data/final_celltype_metadata.txt"
covariate_mat_file="data/cov.mat.txt"
map_covariates_file="data/map_vars_covariates.txt"
expfile="data/57epigenomes.RPKM.pc"
state_annotations_file="data/core_annotation.txt"
generegions_label="gencode_v10"
genefile="data/gencode_genes_full.txt"
set_variables(model, metadatafile,  genefile, covariate_mat_file, map_covariates_file, expfile, state_annotations_file, generegions_label)

## for a particular test
metric="perc"
correction="fdr"
test="wilcox"
property="age"
a.label="Adult"
b.label="Fetal"
pair.label=get_pair_label(a.label, b.label)
celltypes.list=get_celltypes(property, pair.label)
group.a=celltypes.list$a.vec
group.b=celltypes.list$b.vec

#a=get_covariate_mat("sex")
#load("rdata/core_gencode_v10/perc/all_numsonly.rdata")
#fixed.cov.mat=a[colnames(mat),]
#test.mat=mat[1:100,]
#test.totalcounts=totalcounts[1:100]
#test.resids.1=get_residuals_from_mats(fixed.cov.mat, t(test.mat), test.totalcounts, method="log")
#test.resids.2=get_residuals_from_mats(fixed.cov.mat, t(test.mat), test.totalcounts, method="log")
#test.resids.3=get_residuals_from_mats(fixed.cov.mat, t(test.mat), test.totalcounts, method="log")
