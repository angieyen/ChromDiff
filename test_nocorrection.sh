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
for metric in "perc"
do
	for test_type in "wilcox"
	do
		for correction in "fdr"
		do
			for model in "core"
			do
				randomize="FALSE"
				sampled="FALSE"
				expronly="FALSE"
				correctionoff="TRUE"
				a_label="nocorr.${a_option}"
				b_label="nocorr.${b_option}"
				pval_file="all_pvals/${model}/${metric}.${test_type}.${correction}.${property}.${a_label}.${b_label}.txt"
				#if [ ! -e $pval_file ]; then
					bsub -o jobs/nocorr.job -q compbio-week -P compbiofolk -R "rusage[mem=4]" "Rscript test_args.R $property ${a_option} ${b_option} ${a_label} ${b_label} $metric $test_type $correction $model $randomize $sampled $expronly $correctionoff; ./plot_sig_majority_args.sh $property ${a_option} ${b_option} ${a_label} ${b_label} $metric $test_type $correction $model $heightcutoff"
				#fi;
			done
		done
	done
done
