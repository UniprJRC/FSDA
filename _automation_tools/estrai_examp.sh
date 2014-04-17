#!/bin/sh

set -x
#echo "Inserisci il nome dell'archivio FSDA da analizzare : "
#read what

what=`cygpath -a $1`

if [ ! -d "$what" ]
then
        echo "ERRORE :" $what " non esiste !!!!"
        exit 1
fi

mkdir $2/EXAMPLES_test
if [ $? != 0 ]
then
	echo "ERRORE : EXAMPLES_test esiste !!!!"
	exit 1
fi

/usr/bin/find $what -name "*.m" >list_matlab_func

cat list_matlab_func |while read func_file
do
	awk -f /cygdrive/c/Temp/estrai.awk $func_file
done

rm list_matlab_func
echo "I risultati sono in EXAMPLES_test "
exit 0



#echo ""
#echo "Continuo (y/n) ?"
#read risp
#if [ "$risp" != y -a "$risp" != Y ]
#then
#echo " ABORT!!!!! "
#exit 1
#fi
