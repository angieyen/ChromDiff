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
genefile=$1
generegions_label=$2

echo "Processing genefile $genefile..."

mkdir -p genes/${generegions_label}
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5}' $genefile > genes/${generegions_label}/genes.bed
## genenames are sorted lexicographically
awk '(index($1, "#")!=1){print $5}' genes/${generegions_label}/genes.bed | sort -k1,1 > genes/${generegions_label}/genenames.txt
awk '(match($1, "#")==0) {print $6}' $genefile | sort | uniq | wc -l > genes/${generegions_label}/gencode_mapped_symbols.bed
