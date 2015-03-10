#!/bin/sh

set -x
#echo "Inserisci il nome dell'archivio FSDA da analizzare : "
#read what

what=$1

if [ ! -d "$what" ]
then
        echo "ERRORE :" $what " non esiste !!!!"
        exit 1
fi

find $what -name ""*.m"" >lista_func

>func_without_examples.txt
echo -e "The following MATLAB files contain no examples:\n" >>func_without_examples.txt

cat lista_func | while read nome
do
grep "%{" $nome >/dev/null
if [ $? == 1 ]
then
	echo $nome | grep -v "ontents.m" >>func_without_examples.txt
fi
done

echo " Il report si trova in func_without_examples.txt "
unix2dos func_without_examples.txt

rm lista_func
