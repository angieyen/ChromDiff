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
a_label=$2
b_label=$3
test_type=$4
correction=$5
curr_label=$6
heightcutoff=$7
metadatafile=$8
genefile=${9}
covariate_mat_file=${10}
map_covariates_file=${11}
expfile=${12}
state_annotations_file=${13}
generegions_label=${14}

Rscript test_args_v2.R $property ${a_label} ${b_label} $test_type $correction $curr_label $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
sigpval_file="sig_pvals/${curr_label}/perc.${test_type}.${correction}.$property.${a_label}.${b_label}.txt"
if [ -s $sigpval_file ]; then
	## You can first run heightcutoff=0, then pick a height cutoff and re-run plot_sig_majority.sh with the new height cutoff
	./plot_sig_majority_args_v2.sh $property $a_label $b_label $test_type $correction $heightcutoff $curr_label $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
	## calc pvals for expression differences of subgroups
	Rscript calc_pvals_v2.R $property $a_label $b_label $test_type $correction $curr_label $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label


	### use automatically annotated subgroups from mk_clust_heatmaps.R #############
	##### then get MSigDB enrichments
	cd msigdb
	reldir="../"
	Rscript calc_hypergeom_v2.R $property $a_label $b_label $test_type $correction $curr_label ${reldir}$metadatafile ${reldir}$genefile ${reldir}$covariate_mat_file ${reldir}$map_covariates_file ${reldir}$expfile ${reldir}$state_annotations_file $generegions_label

fi
