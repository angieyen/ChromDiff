#!/bin/awk
(NR>2) {len=$3-$2; 
	if ($1=="chrX") {
		x_sum=x_sum+len; 
		x_count=x_count+1;
	} else { 
		sum=sum+len; 
		count=count+1
	} 
} 

END{print "chrX", x_sum, x_count, x_sum/x_count; print "autosomes", sum, count, sum/count;}


