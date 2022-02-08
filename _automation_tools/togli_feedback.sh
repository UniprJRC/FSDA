#!/bin/sh

set -x
echo "Inserisci il nome dell'archivio FSDA da analizzare : "
read what

if [ ! -d "$what" ]
then
        echo "ERRORE :" $what " non esiste !!!!"
        exit 1
fi

qui=`pwd`
mkdir FSDA_senzafeedback
if [ $? != 0 ]
then
	echo "ERRORE : FSDA_senzafeedback esiste !!"
	exit 1
fi

cd $what/helpfiles/FSDA 
find . -name "*.html" >$qui/FSDA_senzafeedback/list_helptomod

cd $qui/FSDA_senzafeedback

cat list_helptomod | while read help_file
do
	awk -f $qui/rm_feedback.awk $what/helpfiles/FSDA/$help_file >$help_file
	unix2dos $help_file

done

rm list_helptomod
echo "I risultati sono in FSDA_senzafeedback "
exit 0



#echo ""
#echo "Continuo (y/n) ?"
#read risp
#if [ "$risp" != y -a "$risp" != Y ]
#then
#echo " ABORT!!!!! "
#exit 1
#fi
