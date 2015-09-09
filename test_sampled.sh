#!/bin/bash
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
property=$1
a_option=$2
b_option=$3
if [ $# -ge 4 ]; then
	heightcutoff=$4
else
	heightcutoff="0"
fi
metric="perc"
test_type="wilcox"
correction="fdr"
statecalls_label="core"
generegions_label="gencode_v10"
curr_label="${statecalls_label}_${generegions_label}"
metadatafile="data/final_celltype_metadata.txt"
genefile="data/gencode_genes_full.txt"
covariate_mat_file="data/cov.mat.txt"
map_covariates_file="data/map_vars_covariates.txt"
expfile="data/57epigenomes.RPKM.pc"
state_annotations_file="data/core_annotation.txt"

numa=8
numb=8
a_label="rand.${a_option}.${numa}"
b_label="rand.${b_option}.${numb}"
pval_file="all_pvals/${model}/${metric}.${test_type}.${correction}.${property}.${a_label}.${b_label}.txt"
if [ ! -e $pval_file ]; then
	bsub -o jobs/rand.job -q compbio-week -P compbiofolk -R "rusage[mem=4]" "./perform_analysis.sh $property $a_label $b_label $test_type $correction $curr_label $heightcutoff $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label " 
fi;
