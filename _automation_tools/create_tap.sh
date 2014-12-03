#!/bin/sh

tapfile=$1/pippo.tap

rm -f $tapfile
awk -f $1/_automation_tools/create_tap.awk execution_log.txt

test_count=`wc -l $tapfile | awk '{ RS=" "; print $1; }'`

echo "1..$test_count" >$1/tempfile.tap
cat $tapfile >>$1/tempfile.tap
rm $tapfile
mv $1/tempfile.tap $tapfile
rm $1/tempfile.tap