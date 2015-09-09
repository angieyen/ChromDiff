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
#for script in "./test_in_parallel.sh" "./test_random.sh" "./test_sampled.sh" "./test_expr.sh" "./test_nocorrection.sh"
for script in "./test_expr.sh" 
do
	#script="./test_in_parallel.sh"
	#script="./test_random.sh"
	#script="./test_sampled.sh"
	#script="./test_expr.sh"


	property="type"
	$script $property "CellLine" "PrimaryCulture" 70
	$script $property "PrimaryCulture" "PrimaryCell" 66
	$script $property "PrimaryCulture" "PrimaryTissue" 66
	$script $property "PrimaryCell" "PrimaryTissue" 83

	property="age"
	$script $property "Adult" "Fetal" 80

	## 3 of the rest 
	property="sex"
	$script $property "Female" "Male" 80 

	property="anatomy"
	$script $property "BRAIN" "MUSCLE" 
	$script $property "BRAIN"  "ESC" 42
	$script $property "BRAIN" "SKIN" 40
	$script $property  "MUSCLE" "ESC" 
	$script $property  "MUSCLE" "SKIN"
	$script $property  "ESC" "SKIN" 40

	property="specialgi"
	$script $property "BRAIN" "GI" 50
	$script $property "SKIN" "GI" 48
	$script $property "ESC" "GI" 45
	$script $property "MUSCLE" "GI"

	######### NEEDS MORE MEMORY ##################
	property="solid_liquid"
	$script $property "SOLID" "LIQUID" 80 
done
	
