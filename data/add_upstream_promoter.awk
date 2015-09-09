BEGIN{
	OFS="\t"
	print "#All gencode protein-coding gene bodies with additional window " window "bp upstream."
	print "#chr", "start", "end", "strand", "gene_id"
	while(getline < "chromInfo.txt") {
		save[$1]=$2
	}
}

(index($1, "#")==0){
	if ($4=="+") {
		start=$2-window
		end=$3
	} else {
		start=$2
		end=$3+window
	}
	
	if (start <0) {
		start=0
	}
	if (end>save[$1]) {
		end=save[$1]
	}
	print $1, start, end, $4, $5, $6
}
