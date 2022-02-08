#!/bin/sh
#
set -x
echo "Devi essere posizionato sopra al folder FSDA"
echo "Continuo ? [y/n]"
read risp
if [ "$risp" == n ]
then
	exit
fi
radice=`pwd`
from=$radice/FSDA
to=$radice/FSDA_modified

mkdir $to
if [ $? != 0 ]
then
	echo "$to ESISTE !!"
	exit 1
fi
workdir="$to/work"
mkdir $workdir

cd $from
find . -name "*.m" -exec grep "docsearch(" {} \; -print |grep -v "docsearch" |grep -v "Contents.m" >$workdir/l_f_docs


cat $workdir/l_f_docs |while read nomefile
do
dname=`dirname "$nomefile"`
mkdir -p $to/$dname
sed -e ' {
	s/docsearch/docsearchFS/
	}' "$nomefile" >$to/"$nomefile"
unix2dos $to/"$nomefile"
done

cd $radice
diff -r FSDA FSDA_modified >$workdir/diff.txt
echo "diff tra FSDA e FSDA_modified in " $workdir/diff.txt

