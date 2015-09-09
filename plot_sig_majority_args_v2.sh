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
## After running mk_clust_heatmaps.R to generate first-time dendrogram/pvclust, you can look for dendrogram with heights in the file pvclust_domstate.pdf in the corresponding plots/core/*/perc/wilcox/fdr/ file
## Then pick a height cutoff and re-run plot_sig_majority.sh with the new height cutoff

test_one=false
if [[ "$test_one" == true ]]; then
	testtype="wilcox"
	correction="fdr"
	property="sex"
	group1="Female"
	group2="Male"
	heightcutoff=83
	model="core_gencode_v10_1kb"
	metadatafile="data/final_celltype_metadata.txt"
	genefile="data/gencode_genes_1kb_up.txt"
	covariate_mat_file="data/cov.mat.txt"
	map_covariates_file="data/map_vars_covariates.txt"
	expfile="data/57epigenomes.RPKM.pc"
	state_annotations_file="data/core_annotation.txt"
	generegions_label="gencode_v10_1kb"
else
	property=$1
	group1=$2
	group2=$3
	testtype=$4
	correction=$5
	heightcutoff=$6
	model=$7
	metadatafile=$8
	genefile=${9}
	covariate_mat_file=${10}
	map_covariates_file=${11}
	expfile=${12}
	state_annotations_file=${13}
	generegions_label=${14}	
fi


label1=$group1
label2=$group2
metric="perc"

sig_dir="sig_pvals"
sig_file="${sig_dir}/${model}/${metric}.${testtype}.${correction}.${property}.${label1}.${label2}.txt"
if [ ! -s $sig_file ]; then
	echo "$sig_file not found or of size 0. Exiting..."
	exit 0
fi
maj_dir="majority"
			
echo "Starting combination of ${model} $metric $testtype $correction $label1 and $label2..."

currprefix="${metric}.${testtype}.${correction}.${property}.${label1}.${label2}"


## what genes had at least one significant feature for this combination?
echo "Identifying genes..."
outdir="sig_majority/${model}"
mkdir -p $outdir
outgenefile="${outdir}/${currprefix}.genes.txt"
awk 'BEGIN{FS="_"} { print $1 }' $sig_file | sort -k1,1 | uniq  > $outgenefile

if [ ! -s $outgenefile ]; then
	echo "No significant genes found for $currprefix." 
	exit 0
else
	subdir="celltype"
	## find majority state for each "significant" gene ###########
	echo "Processing epigenomes..."
	IFS=$'\t'
	header=true
	while read epid statecallfile rest; do
        if $header ;
        then
                header=false
        else
            #echo "Processing $epid features..."
			if [[ "$metric" == "windows" || "$metric" == "deltawindows" ]]; then
				maj_file="${maj_dir}/${model}/${epid}_majstate_windows.txt"
			else
				maj_file="${maj_dir}/${model}/${epid}_majstate.txt"
			fi
			if [ -s $maj_file ]; then
				mkdir -p "${outdir}/${subdir}"
				awk -v outgenefile=${outgenefile} 'BEGIN{ while(getline < outgenefile) { genes[$1]="TRUE"}} 
					{if ($5 in genes){split($6, elts, "_"); found[$5]=elts[1]; delete elts}} 
					END{for (geneid in genes) { if (geneid in found) {print geneid, found[geneid];} else {print geneid, "NA"}}} ' ${maj_file} | sort -k1,1 > "${outdir}/${subdir}/${epid}_${currprefix}.txt"
			fi
        fi
	done < $metadatafile
	
	####### make a single matrix with all the information ########
	echo "Making matrix..." 
	matfile="${outdir}/${currprefix}_matrix.txt"
	## Write column names (genes) 
	awk 'BEGIN{ printf "celltype\t"} {printf "%s\t", $1 } END{print ""}' $outgenefile > $matfile
	## for each celltype and gene, write the majority state
	IFS=$'\t'
	header=true
	while read epid statecallfile rest; do
        if $header ;
        then
                header=false
        else
			currfile="${outdir}/${subdir}/${epid}_${currprefix}.txt"
			if [ -s $currfile ]; then
				awk -v epid=${epid} 'BEGIN{ printf "%s\t", epid} ($2=="NA"){printf "NA\t";} ($2!="NA"){printf "%i\t", $2 } END{print ""}' $currfile >> $matfile
			fi
        fi
	done < $metadatafile

	
	#### with data, generate plots #####
	for plottype in "domstate" 
	do
		echo "Generating plot type $plottype..."
		echo "Generating majority state heatmap..."
		echo Rscript mk_sig_maj_heatmap.R $metric $testtype $correction $property $group1 $group2 $label1 $label2 ${model} $plottype $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label
		Rscript mk_sig_maj_heatmap.R $metric $testtype $correction $property $group1 $group2 $label1 $label2 ${model} $plottype $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label


		echo "Generating expression heatmap..."
		Rscript mk_exp_heatmap.R $metric $testtype $correction $property $group1 $group2 $label1 $label2 ${model} $plottype $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label

		
		echo "Generating clustered heatmap..."
		Rscript mk_clust_heatmaps.R $metric $testtype $correction $property $group1 $group2  $label1 $label2 ${model} $heightcutoff ${plottype} $metadatafile $genefile $covariate_mat_file $map_covariates_file $expfile $state_annotations_file $generegions_label

	done
fi
