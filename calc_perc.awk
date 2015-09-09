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

	ind=1
	genenames_file="genes/" genelabel "/genenames.txt"
	while(getline < genenames_file) {
		genename[ind]=$1
		ind+=1
	}

	currgeneind=1

	countsfile="perc/" model "/numsonly/" celltype "_counts.txt"
}

## moved on to new gene, so print out results for old gene and save new 
($5!=currgeneid) {
	if(currgeneid!="") {	

		## increment through all genenames that had no associated calls (they didn't show up in intersection file)
		## keep printing NAs for these genes
		while(currgeneid!=genename[currgeneind]) {

			## print "NA" for each state that we have
			state=1
			printf "%s\t", genename[currgeneind]
			while(state<=maxstate) {
				printf "NA\t"
				print "0" > countsfile
				state+=1
			}
			print ""
			## increment because we have printed another row
			currgeneind+=1
		}


		currstate=1
		## print each row as geneid covstate1 covstate2 covstate3.... where each column is percent coverage of that gene with the given state
		printf "%s\t", currgeneid
		while(currstate<=maxstate) {
			if(currstate in statecounts) {
				entry=statecounts[currstate]/total
			} else {
				entry=0
			}
			printf "%f\t", entry
			print total > countsfile
			currstate+=1
		}
		print ""	
		currgeneind+=1
	}
	
	## save gene id
	## reset total
	## reset statecounts array
	currgeneid=$5
	total=0
	delete statecounts
}

## save information on this line
($5==currgeneid){
	## calc total length in this gene for this state (geneid in currgeneid and $5, current state in $9)
	split($9, elts, "_")
	state=elts[1]
	newval=($8-$7)
	## take care of cases where ChromHMM region goes beyond gene location
	if ($7<$2) {
		newval-=($2-$7)
	}
	if ($8>$3) {
		newval-=($8-$3)
	}
	statecounts[state]+=newval
	total+=newval
}



END{
	if(currgeneid!="") {	
		## increment through all genenames that had no associated calls (they didn't show up in intersection file)
		## keep printing NAs for these genes
		while(currgeneid!=genename[currgeneind]) {

			## print "NA" for each state that we have
			printf "%s\t", genename[currgeneind]
			state=1
			while(state<=maxstate) {
				printf "NA\t"
				print "0" > countsfile
				state+=1
			}
			print ""
			## increment because we have printed another row
			currgeneind+=1
		}
		
		currstate=1
		## print each row as covstate1 covstate2 covstate3.... where each column is percent coverage of that gene with the given state
		## currgeneid has geneid if necessary: printf "%s\t", currgeneid
		printf "%s\t", currgeneid
		while(currstate<=maxstate) {
			if(currstate in statecounts) {
				entry=statecounts[currstate]/total
			} else {
				entry=0
			}
			printf "%f\t", entry
			print total > countsfile
			currstate+=1
		}
		print ""	
		currgeneind+=1
	}
	
	## we have no more saved information about any genes from input file, so just print "NA" for rest
	## increment through all genenames that had no associated enhancer (they didn't show up in intersection file)
	## keep printing NAs for these genes
	while(currgeneind<=length(genename)) {

	  	## print "NA" for each enhancer state that we have
		printf "%s\t", genename[currgeneind]
		state=1
		while(state<=maxstate) {
			printf "NA\t"
			print "0" > countsfile
			state+=1
		}
		print ""
		## increment because we have printed another row
		currgeneind+=1
	}
}
