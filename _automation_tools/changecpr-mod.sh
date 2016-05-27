#!/bin/sh
#
#set -x

echo "Inserire il path assoluto (Linux) dell'archivio FSDA da modificare : "
echo " es: /home/adminjrc/FSDA_last_checkout "
read from

if [ ! -d "$from" ]
then
  echo "$from non esiste - ERRORE"
  exit 1
fi
echo " Verra' modificata la stringa di Copyright con l'anno corrente e Last Modified con la data di ultima modifica al file .m "
echo " Continuo [s/n] ? \c"
read risp
if [ "$risp" == "" -o "$risp" == "n" ]
then
	echo "Procedura terminata."
	exit 1
fi

to=/tmp/FSDA_modified
mkdir $to
if [ $? != 0 ]
then
  rm -rf /tmp/FSDA_modified
  mkdir $to
fi

cd $from
find . -name "*.m" |grep -v "Contents.m" | grep -v "%" >REP/l_f_copyright

cat REP/l_f_copyright |while read nomefile
do
dname=`dirname "$nomefile"`
mkdir -p $to/$dname

export anno=`date +%Y`

export data=`date -r ${to}/${nomefile}`

vi "$nomefile" <<EOF
:g/% Copyright [0-9,-]*
!!echo $anno
I% Copyright 2008-^[
:g/% Last modified
!!echo $data
I% Last modified ^[:wq
EOF


unix2dos $to/"$nomefile"
done


diff -r $from $to >REP/diff_copyright.txt
echo "Trovi i .m modificati in $to "
echo "Il report della diff tra FSDA e FSDA_modified in " REP/diff_copyright.txt

