#!/bin/awk
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
BEGIN{
	OFS="\t"
	tab="\t"
}

## moved on to new gene, so print out results for old gene and save new 
($5!=last) {
	maxind=0
	max=0
	total=0
	for(ind in statecounts) {
		total+=statecounts[ind]
		if(statecounts[ind]>max) {
			max=statecounts[ind]
			maxind=ind
		}
	}
	if (max!=0) {
		## print gene, majority state, majority state length, total gene length
		printf "%s", geneinfo
		print last, maxind, statecounts[maxind], total
		delete statecounts
	}
	## save gene id and location/strand for later
	last=$5
	geneinfo=$1 tab $2 tab $3 tab $4 tab 
}

## save information on this line
($5==last){
	## calc total length in each gene for each state
	## add chromstate region
	statelength=$8-$7
	statecounts[$9]+=statelength
	## take care of cases where ChromHMM region goes beyond gene location
	##starts before gene, subtract difference
	if ($7<$2) {
		statecounts[$9]-=$2-$7
	}
	## ends after gene, subtract difference
	if ($8>$3) {
		statecounts[$9]-=$8-$3
	}
}
	
## moved on to new gene, so print out results for old gene and save new 
END {
	maxind=0
	max=0
	total=0
	for(ind in statecounts) {
		total+=statecounts[ind]
		if(statecounts[ind]>max) {
			max=statecounts[ind]
			maxind=ind
		}
	}
	if (max!=0) {
		## print gene, majority state, majority state length, total gene length
		printf "%s", geneinfo
		print last, maxind, statecounts[maxind], total
		delete statecounts
	}
}
