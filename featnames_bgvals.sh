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
## check that percentage rows add to 1
#statecalls_label="core"
#generegions_label="gencode_v10"
#states_info="core_annotation.txt"
statecalls_label=$1
generegions_label=$2
states_info=$3

curr_label=${statecalls_label}_${generegions_label}

numstates=`awk '(NR!=1 && $1>max){max=$1} END{print max}' $states_info`
echo "Checking percentage files..."
numerrs_perc=`awk '($1!="NA"){sum=0; i=1; while(i<=NF){sum+=$i; i+=1;}; if (sum<.99999) {print NR, sum}}' perc/${curr_label}/numsonly/*_numsonly.txt | wc -l`
if [ $numerrs_perc -ne 0 ]; then echo "$numerrs_perc errors in $curr_label files" ; exit 1; fi
echo "Writing out featurenames..."
genename_file="genes/${generegions_label}/genenames.txt"
## print all feature names in order from left to right, top to bottom 
mkdir -p featurenames/${curr_label}
awk -f gen_featnames.awk -vmaxstate=${numstates} $genename_file > featurenames/${curr_label}/featurenames.txt
echo "Calculating background values..."
## use all data as background (calc expected (average) background percentage for all gene/chrstate combinations in all celltypes)
mkdir -p backgrounds/${curr_label}
if [ ! -s backgrounds/${curr_label}/perc_background.txt ]; then 
	awk -voutputfile="backgrounds/${curr_label}/perc_background.txt" -f calc_background.awk perc/${curr_label}/numsonly/*_numsonly.txt
fi
## check that background rows add to 1
echo "Checking background files..."
numerrs_bgperc=`awk '($1!="NA"){sum=0; i=1; while(i<=NF){sum+=$i; i+=1;}; if (sum<.99999) {print NR, sum}}' backgrounds/${curr_label}/perc_background.txt | wc -l`
if [ $numerrs_bgperc -ne 0 ]; then echo "$numerrs_bgperc errors in backgrounds/${curr_label}/perc_background.txt"; exit 1; fi

