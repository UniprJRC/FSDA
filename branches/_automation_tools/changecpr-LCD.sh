#!/bin/sh
#
set -x

echo "Inserire il path assoluto (Linux) dell'archivio FSDA da modificare : "
echo " es: /home/adminjrc/FSDA_last_checkout "
read from

if [ ! -d "$from" ]
then
  echo "$from non esiste - ERRORE"
  exit 1
fi
echo " Verra' modificata la stringa di Copyright con l'anno corrente e Last Modified con la LastChangedDate al file .m "
echo " Continuo [s/n] ? \c"
read risp
if [ "$risp" == "" -o "$risp" == "n" ]
then
	echo "Procedura terminata."
	exit 1
fi

currdir=`pwd`

to=/tmp/FSDA_modified
mkdir $to
if [ $? != 0 ]
then
  rm -rf /tmp/FSDA_modified
  mkdir $to
fi

cd $from
find . -name "*.m" |grep -v "Contents.m"  >$currdir/REP/l_f_m

cat $currdir/REP/l_f_m |while read nomefile
do
dname=`dirname "$nomefile"`
mkdir -p $to/$dname


export anno=`date +%Y`

export data=`date -r ${from}/${nomefile} +%d-%m-%Y`

#sed -e 's/% Copyright 2008-.*/% Copyright 2008-'"$anno"'./' -e 's/%Copyright 2008-.*/% Copyright 2008-'"$anno"'./' -e 's/% Last modified .*/%$LastChangedDate::                     #$: Date of the last commit/'  < $from/$nomefile >$to/$nomefile
sed -e 's/% Copyright 2008-.*/% Copyright 2008-'"$anno"'./' -e 's/%Copyright 2008-.*/% Copyright 2008-'"$anno"'./'  < $from/$nomefile >$to/$nomefile
unix2dos $to/"$nomefile"
done


diff -r $from $to >$currdir/REP/diff_cpr_mod.txt
echo "Trovi i .m modificati in $to "
echo "Il report della diff tra FSDA e FSDA_modified in " $currdir/REP/diff_cpr_mod.txt

