#!/bin/sh

#set -x
echo "Inserisci il nome dell'archivio FSDA da analizzare : "
read what

if [ ! -d "$what" ]
then
        echo "ERRORE :" $what " non esiste !!!!"
        exit 1
fi

find $what -name "*.m" >lista_func

>REP/func_without_example.txt

cat lista_func | while read nome
do
grep "%{" $nome >/dev/null
if [ $? == 1 ]
then
	echo $nome | grep -v "ontents.m" >>REP/func_without_example.txt
fi
done

echo " Il report si trova in REP/func_without_example.txt "
unix2dos REP/func_without_example.txt

rm lista_func
