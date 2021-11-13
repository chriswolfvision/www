#!/bin/ksh

if [ $# -ne 2 ]
then
	echo "usage: $0 <filename> <basedir>"
	exit 1
fi
fn=$1
bn=$(basename $fn)
dir=$2

LINES=$(find $dir -name "*.html" -exec grep $bn {} \; | wc -l)

if [ $LINES -le 0 ]
then
	size=$(stat -f %z $fn)
	echo "$size $fn"
fi
