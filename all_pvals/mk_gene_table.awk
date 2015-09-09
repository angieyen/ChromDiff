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
## makes table with one line for each gene seen
BEGIN{OFS="\t"; while(getline < "../genes/gencode_genes_full.bed") {save[$5]=$6; chr[$5]=$1}} 

(NR!=1  && $2<.05){split($1, elts, "_"); id=elts[1]; statenum=elts[2]; 
if(!(id in seen)) {
	print id, statenum, save[id], chr[id], $2;
	seen[id]="Yes"
}
}
