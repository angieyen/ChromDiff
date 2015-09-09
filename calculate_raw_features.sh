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

#statecallfile="statecalls/core/E001_15_coreMarks_mnemonics.bed.gz"
#statecalls_label="core"
#generegions_label="gencode_v10"
#states_info="core_annotation.txt"
epid=$1
statecallfile=$2
statecalls_label=$3
generegions_label=$4
states_info=$5

curr_label=${statecalls_label}_${generegions_label}

## intersect_genes.sh
## overlap gene bodies with chr states
## intersect files are sorted by gene name (5th column) in intersect/${curr_label}/${epid}_genes.txt
## then picks majority state out of considered states and outputs for each gene into majority/${curr_label}/${epid}_majstate.txt 
intfile="intersect/${curr_label}/${epid}.txt"
mkdir -p intersect/${curr_label} majority/${curr_label} 

majfile="majority/${curr_label}/${epid}_majstate.txt"
if [ ! -s $majfile ]; then
	if [[ $statecallfile == *.gz ]]; then
		zcat $statecallfile | intersectBed -wa -wb -a genes/${generegions_label}/genes.bed -b stdin | sort -k5,5 > $intfile ; awk -f pick_majority.awk $intfile > majority/${curr_label}/${epid}_majstate.txt
	else
		intersectBed -wa -wb -a genes/${generegions_label}/genes.bed -b $statecallfile | sort -k5,5 > $intfile ; awk -f pick_majority.awk $intfile > majority/${curr_label}/${epid}_majstate.txt
	fi
fi


## calc_perc_all.sh
## calc value for each gene, scored by % of bases in each state
## outputs in $percfile
percfile="perc/${curr_label}/numsonly/${epid}_numsonly.txt"
if  [ ! -s $percfile ]; then
	mkdir -p perc/${curr_label}/numsonly
	numstates=`awk '(NR!=1 && $1>max){max=$1} END{print max}' $states_info`
	awk -vmaxstate=$numstates -vmodel=${curr_label} -vcelltype=${epid} -vgenelabel=${generegions_label} -f calc_perc.awk $intfile > perc/${curr_label}/${epid}_percentages.txt;
	awk -f get_numsonly.awk perc/${curr_label}/${epid}_percentages.txt > $percfile
	## calc_perc.awk also generates perc/${curr_label}/numsonly/${epid}_counts.txt
fi
