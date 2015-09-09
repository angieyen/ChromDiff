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

for model in "core"
do
	statecallsdir="statecalls/${model}/"
	for int1 in {1..129}
	do
		## make celltype name of the form Exxx as necessary
		if [ $int1 -lt 100 ]; then
    		if [ $int1 -lt 10 ]; 
    			then celltype1="E00${int1}";
    		else 
    			celltype1="E0${int1}";
    		fi
    	else
    		celltype1="E${int1}"
    	fi

		echo "Processing $celltype1... "
		if [ $model == "core" ]; then
			statecallfile="${statecallsdir}${celltype1}_15_coreMarks_mnemonics.bed.gz"
		fi
			
		if [ -s $statecallfile ]
		then
#			if [ ! -s "majority/${model}/${celltype1}_majstate.txt" ]; then
				bodiesintfile="intersect/${model}/${celltype1}_genes.txt"
				## overlap gene bodies with state calls (filter out chrY)
				zcat $statecallfile | awk -f filter_chrY.awk | intersectBed -wa -wb -a genes/gencode_genes.bed -b stdin | sort -k5,5 > $bodiesintfile ; awk -f pick_majority.awk $bodiesintfile > majority/${model}/${celltype1}_majstate.txt 
#			fi
		fi

	done
done
