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
if [ $# -ge 5 ]; then
	memamt=$5
else 
	memamt="4"
fi
randomize="FALSE"
sampled="FALSE"
expronly="FALSE"
nocorrection="FALSE"

for metric in "perc"
do
	for test_type in "wilcox"
	do
		for correction in "fdr" "bonferroni" "BY"
		do
			for model in "core" 
			do
				a_label=$a_option
				b_label=$b_option
				sig_file="sig_pvals/${model}/${metric}.${test_type}.${correction}.${property}.${a_label}.${b_label}.txt"
				if [ ! -e $sig_file ]; then
					currmemamt=10
					bsub -o jobs/test.job -q compbio-week -P compbiofolk -R "rusage[mem=${currmemamt}]"  "Rscript test_args.R $property ${a_option} ${b_option} ${a_label} ${b_label} $metric $test_type $correction $model $randomize $sampled $expronly $nocorrection;  ./plot_sig_majority_args.sh $property $a_option $b_option $a_label ${b_label} $metric $test_type $correction $model $heightcutoff"
				
				else
					currmemamt=$memamt
					#Rscript test_args.R $property ${a_option} ${b_option} ${a_label} ${b_label} $metric $test_type $correction $model $randomize $sampled $expronly $nocorrection;
					bsub -o jobs/test.job -q compbio-week -P compbiofolk -R "rusage[mem=${memamt}]"  "./plot_sig_majority_args.sh $property $a_option $b_option $a_label $b_label $metric $test_type $correction $model $heightcutoff"
					#./plot_sig_majority.sh $sig_file $heightcutoff
 				fi;
			done
		done
	done
done
