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

## must set outputfile (awk -f -voutputfile="background.txt")

{
	## check if this feature was NA for this celltype (if everything was "NA")
	if($1!="NA"){
		## tally # of celltypes that had valid features
		total[FNR]+=1; 
		i=1;
		## sum percent presence of this gene (FNR) in this enhancer (i) 
		while(i<=NF) { save[FNR, i]+=$i; i+=1 }
	}
	if(FNR==1){print "Processing " FILENAME "...";}
}


END{ 
	print "Printing results..."
	line=1;  
	while(line<=FNR) { 
		col=1; 
		while(col<=NF) 
		{
			## if everything was NA, just print all NA's
			if(!(line in total)) {
				printf "NA\t" > outputfile	
			} else {
				# divide sums by counts to obtain average
				printf "%f\t", save[line, col]/total[line] > outputfile
			}
			col+=1;
		} 
		line+=1; 
		print ""> outputfile;
	}
} 
