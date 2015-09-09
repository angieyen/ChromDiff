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
#file=$1
#modeldir=` dirname $file `
#basefile=` basename $file .txt `
#OIFS=$IFS
#IFS="."
#elts=($basefile)
#metric=${elts[0]}
#testtype=${elts[1]}
#correction=${elts[2]}
#property=${elts[3]}
#groupa=${elts[4]}
#groupb=${elts[5]}
#IFS=$OIFS
metric=$1
testtype=$2
correction=$3
property=$4
label1=$5
label2=$6
model=$7

basefile="${metric}.${testtype}.${correction}.${property}.${label1}.${label2}"
infile="../plots/${model}/${label1}.${label2}/${metric}/${testtype}/${correction}/sigfeats_feats_featorder_matched.txt"
echo $infile
if [ -s $infile ]; then	
	echo "Processing $infile ..."
	awk -f mk_gene_table.awk $infile > "${model}/gene_tables/${basefile}_matched.txt"
	echo "${model}/gene_tables/${basefile}_matched.txt"
	awk -f mk_sig_table.awk $infile > "${model}/sig_tables/${basefile}_matched.txt"
fi
