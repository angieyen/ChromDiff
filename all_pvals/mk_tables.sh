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

for model in "core" "core_k27ac" "imputed"
do	
	## only make tables for comparisons with significant results
	sigdir="../sig_pvals/${model}/"
	for file in ${sigdir}*.txt
	do
		basefile=` basename $file .txt `
		OIFS=$IFS
		IFS="."
		elts=($basefile)
		metric=${elts[0]}
		testtype=${elts[1]}
		property=${elts[2]}
		groupa=${elts[3]}
		groupb=${elts[4]}
		IFS=$OIFS
		#infile="../plots/${model}${groupa}.${groupb}/${metric}/${testtype}/sigfeats_feats_featorder_matched.txt"
		infile="${model}/${basefile}.txt"
		if [ -s $sigfile ]; then	
			echo "Processing $infile ..."
			#awk -f mk_table.awk $infile | awk '(! ($1 in seen)){print $0; seen[$1]="Yes"}' > "${model}/tables/${basefile}_matched.txt"
			awk -f mk_gene_table.awk $infile > "${model}/gene_tables/${basefile}_matched.txt"
			awk -f mk_sig_table.awk $infile > "${model}/sig_tables/${basefile}_matched.txt"
		fi
	done
done
