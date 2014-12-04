#!/bin/sh

tapfile=$1/results.tap

rm -f $tapfile

cat execution_log.txt | tr '\\' '/' | sed 's/C://g' | sed 's/://g' >execution_log2.txt

awk -f $1/_automation_tools/create_tap.awk execution_log2.txt

test_count=` grep -n "EXAMPLES_test" execution_log2.txt | wc -l`
# test_count=`wc -l $tapfile | awk '{ RS=" "; print $1; }'`

echo "1..$test_count" >$1/tempfile.tap
cat $tapfile >>$1/tempfile.tap
rm $tapfile
mv $1/tempfile.tap $tapfile
rm -f $1/tempfile.tap