#!/bin/sh

tapfile="pippo.tap"

rm -f $tapfile
awk -f ${WORKSPACE}/_automation_tools/create_tap.awk execution_log.txt

test_count=`wc -l $tapfile | awk '{ RS=" "; print $1; }'`

echo "1..$test_count" >tempfile.tap
cat $tapfile >>tempfile.tap
rm $tapfile
mv tempfile.tap $tapfile
rm tempfile.tap