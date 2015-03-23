#!/bin/sh

set -x
>REP/help_href-references.txt
>REP/help_href-ERROR.txt

echo "Inserisci il nome dell'archivio FSDA da analizzare : "
read what

if [ ! -d "$what" ]
then
        echo "ERRORE :" $what " non esiste !!!!"
        exit 1
fi

mydir=`pwd`
cd $what/helpfiles/FSDA

find . -name "*.html" |while read help_file
do
	awk -f $mydir/check_hr.awk $help_file >>$mydir/REP/help_href-references.txt
done


cat $mydir/REP/help_href-references.txt |while read href
do
 	init=${href:0:1}
	if [ $init ==  "." ]
	then 
		dove=$href
	fi
	if [ ! -e "$href" ]
	then	
		echo "ERRORE in $dove : href=" $href " --- IL FILE NON ESISTE" >>$mydir/REP/help_href-ERROR.txt
	fi
done

unix2dos $mydir/REP/help_href-references.txt
unix2dos $mydir/REP/help_href-ERROR.txt

echo "I risultati sono in REP/help_href-references.txt e REP/help_href-ERROR.txt"
exit 0



#echo ""
#echo "Continuo (y/n) ?"
#read risp
#if [ "$risp" != y -a "$risp" != Y ]
#then
#echo " ABORT!!!!! "
#exit 1
#fi
