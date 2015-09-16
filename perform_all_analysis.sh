#!/bin/bash

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


statecalls_label="core"
generegions_label="gencode_v10"
curr_label="${statecalls_label}_${generegions_label}"
heightcutoff=0
metadatafile="data/final_celltype_metadata.txt"
genefile="data/gencode_genes_full.txt"
covariate_mat_file="data/cov.mat.txt"
map_covariates_file="data/map_vars_covariates.txt"
expfile="data/57epigenomes.RPKM.pc"
state_annotations_file="data/core_annotation.txt"

correction="fdr"
qsub -cwd -j y -b y -q long -o jobs/CellLine_PrimaryCulture.job ./perform_analysis.sh "type" "CellLine" "PrimaryCulture" "wilcox" $correction $curr_label 39 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/PrimaryCulture_PrimaryCell.job ./perform_analysis.sh "type" "PrimaryCulture" "PrimaryCell" "wilcox" $correction $curr_label 69 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/PrimaryCulture_PrimaryTissue.job ./perform_analysis.sh "type" "PrimaryCulture" "PrimaryTissue" "wilcox" $correction $curr_label 83 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/PrimaryCell_PrimaryTissue.job ./perform_analysis.sh "type" "PrimaryCell" "PrimaryTissue" "wilcox" $correction $curr_label 83 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label

qsub -cwd -j y -b y -q long -o jobs/Adult_Fetal.job ./perform_analysis.sh "age" "Adult" "Fetal" "wilcox" $correction $curr_label 80 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label

qsub -cwd -j y -b y -q long -o jobs/Female_Male.job ./perform_analysis.sh "sex" "Female" "Male" "wilcox" $correction $curr_label 80 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label

qsub -cwd -j y -b y -q long -o jobs/BRAIN_MUSCLE.job ./perform_analysis.sh "anatomy" "BRAIN" "MUSCLE" "wilcox" $correction $curr_label 0 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/BRAIN_ESC.job ./perform_analysis.sh "anatomy" "BRAIN" "ESC" "wilcox" $correction $curr_label 49 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/BRAIN_SKIN.job ./perform_analysis.sh "anatomy" "BRAIN" "SKIN" "wilcox" $correction $curr_label 44 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/MUSCLE_ESC.job ./perform_analysis.sh "anatomy" "MUSCLE" "ESC" "wilcox" $correction $curr_label 0 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/MUSCLE_SKIN.job ./perform_analysis.sh "anatomy" "MUSCLE" "SKIN" "wilcox" $correction $curr_label 0 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/ESC_SKIN.job ./perform_analysis.sh "anatomy" "ESC" "SKIN" "wilcox" $correction $curr_label 0 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label


qsub -cwd -j y -b y -q long -o jobs/BRAIN_GI.job ./perform_analysis.sh "specialgi" "BRAIN" "GI" "wilcox" $correction $curr_label 50 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/SKIN_GI.job ./perform_analysis.sh "specialgi" "SKIN" "GI" "wilcox" $correction $curr_label 48 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/ESC_GI.job ./perform_analysis.sh "specialgi" "ESC" "GI" "wilcox" $correction $curr_label 48 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
qsub -cwd -j y -b y -q long -o jobs/MUSCLE_GI.job ./perform_analysis.sh "specialgi" "MUSCLE" "GI" "wilcox" $correction $curr_label 0 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label

qsub -cwd -j y -b y -q long -o jobs/SOLID_LIQUID.job ./perform_analysis.sh "solid_liquid" "SOLID" "LIQUID" "wilcox" $correction $curr_label 88 $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
