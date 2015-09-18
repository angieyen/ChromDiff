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

## SEE README.TXT FOR FILE FORMAT DETAILS 

### PART 1: PREPROCESSING ###
### VARIABLES, LABELS, AND PARAMETERS TO SET ###
metadatafile="data/final_celltype_metadata.txt"
genefile="data/gencode_genes_full.txt"
generegions_label="gencode_v10"
statecalls_label="core"
states_info="data/core_annotation.txt"
covariate_mat_file="data/cov.mat.txt"
map_covariates_file="data/map_vars_covariates.txt"
expfile="data/57epigenomes.RPKM.pc"
state_annotations_file="data/core_annotation.txt"
### END OF VARIABLES TO SET ####

## Step 1: process and format gene info
./process_gene_info.sh $genefile $generegions_label
 
## Step 2: Calculate feature values
./calculate_all_raw_features.sh $metadatafile $statecalls_label $generegions_label $states_info
## NOTE: can also use ./calculate_raw_features.sh on each epid separately to parallelize


## Step 3: check feature values, write out feature names, calculate background vals, check background values
./featnames_bgvals.sh $statecalls_label $generegions_label $states_info


## PART 2: COMPARING EPIGENOMIC GROUPS ###
## COMPARATIVE GROUPS AND PARAMETERS TO SET ###
property="sex"
a_option="Female"
b_option="Male"
test_type="wilcox"
correction="fdr"
heightcutoff=0
### END OF VARIABLES TO SET ####

## Step 4: Run ChromDiff on a particular epigenomic group comparison
## Note: can run analysis the first time with default height cutoff of 0 for clustering 
## then look at plots/${curr_label}/${a_option}.${b_option}/perc/${test_type}/${correction}/hclust_domstate.pdf to choose manual cutoff
curr_label="${statecalls_label}_${generegions_label}"
./perform_analysis.sh $property $a_option $b_option $test_type $correction $curr_label $heightcutoff $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label 
