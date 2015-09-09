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
	if [ $model == "core" ]; 
		then numstates=15
	fi
	
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

		inputfile1=intersect/${model}/${celltype1}_genes.txt
		#Calc percentages for gene bodies
		if [ -s ${inputfile1} ]; then
    		awk -vmaxstate=$numstates -vmodel=${model} -vcelltype=${celltype1} -f calc_perc.awk ${inputfile1} > perc/${model}/${celltype1}_percentages.txt; 
    		awk -f get_numsonly.awk perc/${model}/${celltype1}_percentages.txt > perc/${model}/numsonly/${celltype1}_numsonly.txt
			#mv perc/${model}/${celltype1}_counts.txt perc/${model}/numsonly/${celltype1}_counts.txt
		fi
	done

done
