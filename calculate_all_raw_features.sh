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
metadatafile=$1
statecalls_label=$2
generegions_label=$3
states_info=$4

IFS=$'\t'
header=true
while read epid statecallfile rest; do
	if $header ;
	then
		header=false
	else
		echo "Processing $epid features..."
		./calculate_raw_features.sh $epid $statecallfile $statecalls_label $generegions_label $states_info
	fi
done < $metadatafile
