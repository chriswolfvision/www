#!/bin/ksh

if [ $# -ne 1 ]
then
	echo "usage: $0 <directory>"
	echo "(list of files to check read from stdin)"
	exit 1
fi
dir=$1

while read i
do 	
	# echo "checking for $fname"
	./scripts/isrefed.sh $i $dir
done
